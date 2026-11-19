/**
 * Extinct Chambered Cephalopod Hydrostatic Simulation
 * 
 * Biomechanical & Paleobiological Modeling:
 * - Trueman, A. E. (1941) Q. J. Geol. Soc. Lond. 96: 339-383.
 * - Raup, D. M. (1966) J. Paleontol. 40: 1178-1190.
 * - Westermann, G. E. G. (1977) Paleobiology 3: 300-321.
 * - Saunders, W. B. & Shapiro, E. A. (1986) Paleobiology 12: 64-79.
 * - Peterman, D. J. et al. (2020) Sci. Rep. 10: 1138.
 */

// Global state and simulation parameters
const SIM_PARAMS = {
  // Shell morphology
  shellLength: 600,            // L in mm
  curvature: 0.0,              // C in [0, 1] (0 = orthocone, 1 = nautilocone)
  apertureDiam: 80,            // D_aperture in mm
  apexDiam: 4,                 // D_apex in mm
  livingChamberLength: 180,    // L_bc in mm
  
  // Shell wall & septa
  shellThickness: 2.0,         // t_shell in mm
  septaThickness: 1.0,         // t_septum in mm
  chamberDist: 10,             // delta_s in mm

  // Siphuncular anatomy
  siphuncleDiam: 10,           // D_siph in mm (diameter at living chamber septum)
  siphunclePos: 50,            // % position: 0% = Ventral margin, 50% = Center, 100% = Dorsal margin
  siphuncleDensity: 2.70,      // rho_siph in g/cm3 (defaults to shell density)
  syncSiphDensity: true,       // whether siphuncle density tracks shell density

  // Hyponome / siphon jet (jet propulsion)
  siphonAngle: 90,             // degrees: angle of siphon jet relative to aperture axis (0° = forward, 90° = perpendicular, 180° = backward)
  siphonForce: 1.0,            // N: magnitude of jet propulsion impulse

  // Densities & Ballast
  shellDensity: 2.70,          // rho_shell in g/cm3 (calcite=2.71, aragonite=2.93)
  chamberGasDensity: 0.0012,   // rho_chamber in g/cm3 (Air=0.0012, CO2=0.00198, Water=1.000)
  weightFactor: 1.00,          // factor for soft tissue (~1.05 g/cm3)
  cameralFluid: 0,             // % residual fluid in camerae
  waterDensity: 1.025,         // rho_water in g/cm3 (marine=1.025, fresh=1.000)

  // Visual options
  showCoM: true,
  showCoB: true,
  showMoI: true,
  showForces: true,
  waterDamping: true
};

// Physics constants
const GRAVITY = 9.80665; // m/s^2

/**
 * -------------------------------------------------------------
 * 1. Conch Geometric & Mass Properties Generator
 * -------------------------------------------------------------
 */
class CephalopodModel {
  constructor(params) {
    this.update(params);
  }

  update(params) {
    this.p = Object.assign({}, params);
    this.L = this.p.shellLength; // mm
    this.C = this.p.curvature;
    this.D_ap = this.p.apertureDiam;
    this.D_apex = this.p.apexDiam;
    this.L_bc = Math.min(this.p.livingChamberLength, this.L * 0.95);
    this.t_shell = Math.min(this.p.shellThickness, this.D_ap * 0.4);
    this.t_septum = this.p.septaThickness;
    this.delta_s = Math.max(2, this.p.chamberDist);
    this.rho_shell = this.p.shellDensity; // g/cm3
    this.siphDiam = Math.max(0.5, this.p.siphuncleDiam !== undefined ? this.p.siphuncleDiam : 10);
    this.siphPos = Math.max(0, Math.min(100, this.p.siphunclePos !== undefined ? this.p.siphunclePos : 50)) / 100; // 0 = ventral, 0.5 = center, 1.0 = dorsal
    this.siphDensity = this.p.siphuncleDensity !== undefined ? this.p.siphuncleDensity : this.rho_shell;
    this.rho_chamber = (this.p.chamberGasDensity !== undefined) ? this.p.chamberGasDensity : 0.0012; // g/cm3 (Air, CO2, or liquid)
    this.tissueDensity = 1.05 * this.p.weightFactor; // g/cm3
    this.rho_water = this.p.waterDensity;
    this.cameralFluidFrac = this.p.cameralFluid / 100;

    this.numSlices = 300;
    this.ds = this.L / this.numSlices;

    this.buildGeometry();
    this.computeMassProperties();
  }

  buildGeometry() {
    this.slices = [];
    this.centerline = [];
    this.ventralWall = [];
    this.dorsalWall = [];
    this.innerVentral = [];
    this.innerDorsal = [];

    // Continuous logarithmic curvature formulation
    // When C = 0: straight cone (orthocone)
    // When C > 0: tangent angle phi(s) grows logarithmically, coiled planispiral
    const maxTheta = 3.6 * Math.PI; // ~1.8 whorls
    const beta = 2.4;
    const totalPhi = this.C * maxTheta;

    let curX = 0;
    let curY = 0;

    for (let i = 0; i <= this.numSlices; i++) {
      const s = (i / this.numSlices) * this.L; // arc length from apex [0, L]
      const u = s / this.L;

      // Tangent angle
      let phi = 0;
      if (this.C > 0) {
        phi = totalPhi * Math.log(1 + beta * u) / Math.log(1 + beta);
      }

      const cosP = Math.cos(phi);
      const sinP = Math.sin(phi);

      // Normal vector pointing to ventral side
      const nx = -sinP;
      const ny = cosP;

      // Outer radius and inner radius (conical taper)
      const r_out = (this.D_apex / 2) + ((this.D_ap / 2) - (this.D_apex / 2)) * u;
      const r_in = Math.max(0.1, r_out - this.t_shell);

      if (i > 0) {
        curX += cosP * this.ds;
        curY += sinP * this.ds;
      }

      const ptCenter = { x: curX, y: curY, s, phi, r_out, r_in, nx, ny, cosP, sinP };
      this.centerline.push(ptCenter);

      // Boundary outlines
      this.ventralWall.push({ x: curX + r_out * nx, y: curY + r_out * ny });
      this.dorsalWall.push({ x: curX - r_out * nx, y: curY - r_out * ny });
      this.innerVentral.push({ x: curX + r_in * nx, y: curY + r_in * ny });
      this.innerDorsal.push({ x: curX - r_in * nx, y: curY - r_in * ny });

      if (i < this.numSlices) {
        const sMid = s + this.ds * 0.5;
        const uMid = sMid / this.L;
        const r_out_mid = (this.D_apex / 2) + ((this.D_ap / 2) - (this.D_apex / 2)) * uMid;
        const r_in_mid = Math.max(0.1, r_out_mid - this.t_shell);
        const isLivingChamber = sMid >= (this.L - this.L_bc);

        const sliceData = {
          index: i,
          s: sMid,
          u: uMid,
          r_out: r_out_mid,
          r_in: r_in_mid,
          isLivingChamber: isLivingChamber,
          x: curX + cosP * (this.ds * 0.5),
          y: curY + sinP * (this.ds * 0.5),
          phi: phi,
          nx: nx,
          ny: ny
        };

        // If in phragmocone, compute siphuncle cone radius and position
        if (!isLivingChamber) {
          const L_phrag = Math.max(1, this.L - this.L_bc);
          const u_phrag = Math.min(1.0, sMid / L_phrag);
          const r_siph_base = this.siphDiam * 0.5;
          const r_siph_apex = Math.max(0.2, r_siph_base * (this.D_apex / this.D_ap));
          const r_siph = Math.min(r_in_mid - 0.2, Math.max(0.1, r_siph_apex + (r_siph_base - r_siph_apex) * u_phrag));
          const r_clear = Math.max(0, r_in_mid - r_siph);
          // 0% = Ventral (+nx), 50% = Center (0), 100% = Dorsal (-nx)
          const offset = r_clear * (1 - 2 * this.siphPos);

          sliceData.r_siph = r_siph;
          sliceData.siph_x = sliceData.x + offset * nx;
          sliceData.siph_y = sliceData.y + offset * ny;
        }

        this.slices.push(sliceData);
      }
    }

    // Siphuncular cone boundary curves for 2D & 3D rendering
    this.siphuncle = [];
    const L_phrag = Math.max(1, this.L - this.L_bc);
    const r_siph_base = this.siphDiam * 0.5;
    const r_siph_apex = Math.max(0.2, r_siph_base * (this.D_apex / this.D_ap));

    for (let i = 0; i <= this.numSlices; i++) {
      const pt = this.centerline[i];
      if (pt.s > (this.L - this.L_bc)) break;

      const u_phrag = Math.min(1.0, pt.s / L_phrag);
      const r_siph = Math.min(pt.r_in - 0.2, Math.max(0.1, r_siph_apex + (r_siph_base - r_siph_apex) * u_phrag));
      const r_clear = Math.max(0, pt.r_in - r_siph);
      const offset = r_clear * (1 - 2 * this.siphPos);

      const cx = pt.x + offset * pt.nx;
      const cy = pt.y + offset * pt.ny;

      this.siphuncle.push({
        s: pt.s,
        cx,
        cy,
        r_siph,
        vx: cx + r_siph * pt.nx,
        vy: cy + r_siph * pt.ny,
        dx: cx - r_siph * pt.nx,
        dy: cy - r_siph * pt.ny
      });
    }

    // Septa positions along phragmocone (s < L - L_bc)
    this.septa = [];
    let sepS = this.L - this.L_bc;
    while (sepS > this.D_apex) {
      // Find interpolated position along centerline
      const sliceIdx = Math.min(this.numSlices - 1, Math.floor((sepS / this.L) * this.numSlices));
      const pt = this.centerline[sliceIdx];
      const r_in = Math.max(0.1, pt.r_out - this.t_shell);

      const u_phrag = Math.min(1.0, sepS / L_phrag);
      const r_siph = Math.min(r_in - 0.2, Math.max(0.1, r_siph_apex + (r_siph_base - r_siph_apex) * u_phrag));

      this.septa.push({
        s: sepS,
        x: pt.x,
        y: pt.y,
        r_in: r_in,
        r_siph: r_siph,
        nx: pt.nx,
        ny: pt.ny,
        phi: pt.phi
      });

      sepS -= this.delta_s;
    }
  }

  computeMassProperties() {
    let totalMass = 0; // grams
    let sumMX = 0;
    let sumMY = 0;

    let shellMass = 0;
    let septaMass = 0;
    let bodyMass = 0;
    let fluidMass = 0;
    let siphMass = 0;
    let totalOuterVol = 0; // cm3

    // 1. Conch wall, living chamber, cameral matter, and siphuncle slices
    for (const slice of this.slices) {
      // Volumes in mm3: dV = pi * (r_out^2 - r_in^2) * ds
      const dV_outer_mm3 = Math.PI * (slice.r_out * slice.r_out) * this.ds;
      const dV_inner_mm3 = Math.PI * (slice.r_in * slice.r_in) * this.ds;
      const dV_wall_mm3 = dV_outer_mm3 - dV_inner_mm3;

      totalOuterVol += dV_outer_mm3 * 1e-3; // cm3

      // Conch wall mass
      const dm_wall = (dV_wall_mm3 * 1e-3) * this.rho_shell; // grams
      shellMass += dm_wall;
      totalMass += dm_wall;
      sumMX += slice.x * dm_wall;
      sumMY += slice.y * dm_wall;

      if (slice.isLivingChamber) {
        // Soft visceral body tissue in living chamber
        const dm_body = (dV_inner_mm3 * 1e-3) * this.tissueDensity; // grams
        bodyMass += dm_body;
        totalMass += dm_body;
        sumMX += slice.x * dm_body;
        sumMY += slice.y * dm_body;
      } else {
        // Siphuncle conical matter in this slice
        const dV_siph_mm3 = Math.PI * (slice.r_siph * slice.r_siph) * this.ds;
        const dm_siph = (dV_siph_mm3 * 1e-3) * this.siphDensity;
        siphMass += dm_siph;
        totalMass += dm_siph;
        sumMX += slice.siph_x * dm_siph;
        sumMY += slice.siph_y * dm_siph;

        // Remaining cameral matter (gas / liquid)
        const dV_cameral_mm3 = Math.max(0, dV_inner_mm3 - dV_siph_mm3);
        const effectiveChamberDensity = this.cameralFluidFrac * this.rho_water + (1 - this.cameralFluidFrac) * this.rho_chamber;
        const dm_cameral = (dV_cameral_mm3 * 1e-3) * effectiveChamberDensity;
        fluidMass += dm_cameral;
        totalMass += dm_cameral;
        sumMX += slice.x * dm_cameral;
        sumMY += slice.y * dm_cameral;
      }
    }

    // 2. Septa (with siphuncular foramen perforation subtracted)
    for (const sep of this.septa) {
      const A_sep_mm2 = Math.max(0, Math.PI * (sep.r_in * sep.r_in - sep.r_siph * sep.r_siph));
      const V_sep_mm3 = A_sep_mm2 * this.t_septum;
      const dm_sep = (V_sep_mm3 * 1e-3) * this.rho_shell;
      septaMass += dm_sep;
      totalMass += dm_sep;
      sumMX += sep.x * dm_sep;
      sumMY += sep.y * dm_sep;
    }

    this.totalMass = totalMass; // grams
    this.shellMass = shellMass;
    this.septaMass = septaMass;
    this.bodyMass = bodyMass;
    this.fluidMass = fluidMass;
    this.siphMass = siphMass;
    this.totalOuterVolume = totalOuterVol; // cm3

    // Center of Mass (CoM / CG) in local coordinates (mm)
    this.com = {
      x: sumMX / totalMass,
      y: sumMY / totalMass
    };

    // 3. 2D Moment of Inertia Tensor relative to CoM
    let Ixx = 0; // g * mm2
    let Iyy = 0;
    let Ixy = 0;

    for (const slice of this.slices) {
      const dx_w = slice.x - this.com.x;
      const dy_w = slice.y - this.com.y;

      const dV_outer_mm3 = Math.PI * (slice.r_out * slice.r_out) * this.ds;
      const dV_inner_mm3 = Math.PI * (slice.r_in * slice.r_in) * this.ds;
      const dm_wall = (dV_outer_mm3 - dV_inner_mm3) * 1e-3 * this.rho_shell;

      Ixx += dm_wall * (dy_w * dy_w);
      Iyy += dm_wall * (dx_w * dx_w);
      Ixy += dm_wall * (dx_w * dy_w);

      if (slice.isLivingChamber) {
        const dm_body = (dV_inner_mm3 * 1e-3) * this.tissueDensity;
        Ixx += dm_body * (dy_w * dy_w);
        Iyy += dm_body * (dx_w * dx_w);
        Ixy += dm_body * (dx_w * dy_w);
      } else {
        // Siphuncle inertia
        const dV_siph_mm3 = Math.PI * (slice.r_siph * slice.r_siph) * this.ds;
        const dm_siph = (dV_siph_mm3 * 1e-3) * this.siphDensity;
        const dx_s = slice.siph_x - this.com.x;
        const dy_s = slice.siph_y - this.com.y;
        Ixx += dm_siph * (dy_s * dy_s);
        Iyy += dm_siph * (dx_s * dx_s);
        Ixy += dm_siph * (dx_s * dy_s);

        // Cameral matter inertia
        const dV_cameral_mm3 = Math.max(0, dV_inner_mm3 - dV_siph_mm3);
        const effectiveChamberDensity = this.cameralFluidFrac * this.rho_water + (1 - this.cameralFluidFrac) * this.rho_chamber;
        const dm_cameral = (dV_cameral_mm3 * 1e-3) * effectiveChamberDensity;
        Ixx += dm_cameral * (dy_w * dy_w);
        Iyy += dm_cameral * (dx_w * dx_w);
        Ixy += dm_cameral * (dx_w * dy_w);
      }
    }

    for (const sep of this.septa) {
      const dx = sep.x - this.com.x;
      const dy = sep.y - this.com.y;
      const A_sep_mm2 = Math.max(0, Math.PI * (sep.r_in * sep.r_in - sep.r_siph * sep.r_siph));
      const V_sep_mm3 = A_sep_mm2 * this.t_septum;
      const dm_sep = (V_sep_mm3 * 1e-3) * this.rho_shell;

      Ixx += dm_sep * (dy * dy);
      Iyy += dm_sep * (dx * dx);
      Ixy += dm_sep * (dx * dy);
    }

    this.Izz = Ixx + Iyy; // Polar moment of inertia (g * mm2)
    // Convert to SI units (kg * m2): 1 g*mm2 = 1e-3 kg * 1e-6 m2 = 1e-9 kg*m2
    this.Izz_SI = this.Izz * 1e-9;

    // Principal Moments of Inertia and Principal Axes
    const trace = Ixx + Iyy;
    const diff = Ixx - Iyy;
    const disc = Math.sqrt(diff * diff + 4 * Ixy * Ixy);
    const lambda1 = (trace + disc) / 2;
    const lambda2 = (trace - disc) / 2;

    this.principalAngle = 0.5 * Math.atan2(2 * Ixy, diff);
    this.k1 = Math.sqrt(lambda1 / totalMass); // Principal radius of gyration (mm)
    this.k2 = Math.sqrt(lambda2 / totalMass); // Secondary radius of gyration (mm)
    this.k_gyr = Math.sqrt(this.Izz / totalMass); // Polar radius of gyration (mm)
  }
}

/**
 * -------------------------------------------------------------
 * 2. Hydrostatic Solver & Equilibrium Calculator
 * -------------------------------------------------------------
 */
class HydrostaticSolver {
  constructor(model) {
    this.model = model;
  }

  /**
   * Calculates submerged volume, buoyant force, center of buoyancy (CoB),
   * and restoring torque for a given shell pose in world coordinates.
   * Water surface is at Y = 0.
   */
  evaluatePose(posX, posY, theta) {
    const cosT = Math.cos(theta);
    const sinT = Math.sin(theta);
    const com = this.model.com;

    let subVol_cm3 = 0;
    let sumBx = 0;
    let sumBy = 0;

    for (const slice of this.model.slices) {
      // Local slice position relative to CoM
      const lx = slice.x - com.x;
      const ly = slice.y - com.y;

      // Slice center in world coordinates
      const wx = posX + (lx * cosT - ly * sinT);
      const wy = posY + (lx * sinT + ly * cosT);

      const R = slice.r_out;
      const depth = -wy; // depth below water line Y=0

      let f_sub = 0;
      let y_offset = 0;

      if (depth >= R) {
        // Entirely submerged
        f_sub = 1.0;
        y_offset = 0;
      } else if (depth <= -R) {
        // Entirely above water
        f_sub = 0.0;
      } else {
        // Partially submerged circular cross-section
        const h = depth + R; // submerged depth from bottom
        const h_norm = Math.max(0, Math.min(2 * R, h));
        const theta_sub = 2 * Math.acos(Math.max(-1, Math.min(1, 1 - h_norm / R)));
        const area_sub = 0.5 * R * R * (theta_sub - Math.sin(theta_sub));
        f_sub = area_sub / (Math.PI * R * R);

        // Centroid of circular segment
        const sinHalf = Math.sin(theta_sub / 2);
        const denom = (theta_sub - Math.sin(theta_sub));
        if (denom > 1e-6) {
          const c_from_center = (4 * R * sinHalf * sinHalf * sinHalf) / (3 * denom);
          y_offset = -c_from_center;
        }
      }

      if (f_sub > 0) {
        const dV_cm3 = (Math.PI * R * R * this.model.ds * 1e-3) * f_sub;
        const bx = wx;
        const by = wy + y_offset;

        subVol_cm3 += dV_cm3;
        sumBx += bx * dV_cm3;
        sumBy += by * dV_cm3;
      }
    }

    const buoyantForce = (subVol_cm3 * 1e-3) * (this.model.rho_water * 1000) * (GRAVITY * 1e-3); // Newtons
    const gravityForce = (this.model.totalMass * 1e-3) * GRAVITY; // Newtons

    let cobX = posX;
    let cobY = posY;

    if (subVol_cm3 > 1e-6) {
      cobX = sumBx / subVol_cm3;
      cobY = sumBy / subVol_cm3;
    }

    // Restoring torque about CoM: tau = (r_B - r_G) x F_B
    // In 2D: tau = (cobX - posX) * F_B (in N * mm) -> N * m
    const armX_m = (cobX - posX) * 1e-3;
    const torque = armX_m * buoyantForce; // N * m

    return {
      subVol_cm3,
      buoyantForce,
      gravityForce,
      cob: { x: cobX, y: cobY },
      torque,
      netForceY: buoyantForce - gravityForce
    };
  }

  /**
   * Find vertical draft Y where buoyantForce == gravityForce for angle theta
   */
  findEquilibriumDraft(theta) {
    let yLow = -this.model.L * 1.5;
    let yHigh = this.model.L * 1.5;

    for (let iter = 0; iter < 28; iter++) {
      const yMid = (yLow + yHigh) * 0.5;
      const res = this.evaluatePose(0, yMid, theta);

      if (res.netForceY < 0) {
        // Sinks too much -> need to submerge more (lower Y)
        yHigh = yMid;
      } else {
        yLow = yMid;
      }
    }

    return (yLow + yHigh) * 0.5;
  }

  /**
   * Solves for the static equilibrium floating orientation and stability index
   */
  solveEquilibrium() {
    const maxDisplacementMass = this.model.totalOuterVolume * this.model.rho_water; // g
    const isFloating = maxDisplacementMass >= this.model.totalMass;

    if (!isFloating) {
      return {
        isFloating: false,
        angleDeg: 0,
        draft: -100,
        com: { x: 0, y: -100 },
        cob: { x: 0, y: -100 },
        stabilityIndex: 0,
        maxDisplacementMass
      };
    }

    // Sample torque across 360 degrees to locate stable zero crossing
    const numSamples = 72;
    let bestAngle = 0;
    let minTorqueAbs = Infinity;

    // Scan for zero crossing of restoring torque
    let prevTorque = null;
    let prevAngle = 0;
    let candidateAngles = [];

    for (let i = 0; i <= numSamples; i++) {
      const angle = -Math.PI + (i / numSamples) * 2 * Math.PI;
      const draft = this.findEquilibriumDraft(angle);
      const res = this.evaluatePose(0, draft, angle);

      if (prevTorque !== null) {
        // Check for zero crossing with negative slope (stable restoring equilibrium)
        if (prevTorque > 0 && res.torque <= 0) {
          candidateAngles.push((prevAngle + angle) * 0.5);
        }
      }

      if (Math.abs(res.torque) < minTorqueAbs) {
        minTorqueAbs = Math.abs(res.torque);
        bestAngle = angle;
      }

      prevTorque = res.torque;
      prevAngle = angle;
    }

    let eqAngle = candidateAngles.length > 0 ? candidateAngles[0] : bestAngle;

    // Refine eqAngle with bisection
    let aLow = eqAngle - 0.15;
    let aHigh = eqAngle + 0.15;

    for (let iter = 0; iter < 16; iter++) {
      const aMid = (aLow + aHigh) * 0.5;
      const draft = this.findEquilibriumDraft(aMid);
      const res = this.evaluatePose(0, draft, aMid);

      if (res.torque > 0) {
        aLow = aMid;
      } else {
        aHigh = aMid;
      }
    }

    eqAngle = (aLow + aHigh) * 0.5;
    const eqDraft = this.findEquilibriumDraft(eqAngle);
    const finalPose = this.evaluatePose(0, eqDraft, eqAngle);

    // Saunders & Shapiro (1986) Stability Index: St = |r_B - r_G| / D_shell
    const sepDistance = Math.hypot(finalPose.cob.x, finalPose.cob.y - eqDraft);
    const stabilityIndex = sepDistance / this.model.D_ap;

    return {
      isFloating: true,
      angleDeg: (eqAngle * 180 / Math.PI),
      angleRad: eqAngle,
      draft: eqDraft,
      com: { x: 0, y: eqDraft },
      cob: finalPose.cob,
      sepDistance,
      stabilityIndex,
      buoyantForce: finalPose.buoyantForce,
      gravityForce: finalPose.gravityForce,
      maxDisplacementMass
    };
  }
}

/**
 * -------------------------------------------------------------
 * 3. Dynamic Physics & Interactive Spring Drag Engine
 * -------------------------------------------------------------
 */
class PhysicsEngine {
  constructor(model, hydro) {
    this.model = model;
    this.hydro = hydro;
    this.resetToEquilibrium();

    // Mouse drag interaction state
    this.isDragging = false;
    this.dragPickBody = { x: 0, y: 0 }; // picked point relative to CoM (mm)
    this.mouseWorld = { x: 0, y: 0 };
    this.appliedForce = { x: 0, y: 0, mag: 0 };
    this.appliedTorque = 0;

    // Siphon jet impulse state
    this.siphonImpulse = { x: 0, y: 0, mag: 0, ttl: 0 }; // ttl = time-to-live in seconds
    this.siphonRefAngleBody = undefined; // computed on first step() call
  }

  resetToEquilibrium() {
    this.eq = this.hydro.solveEquilibrium();
    this.posX = 0;
    this.posY = this.eq.draft;
    this.angle = this.eq.isFloating ? this.eq.angleRad : 0;
    this.vx = 0;
    this.vy = 0;
    this.omega = 0;
    this.isSettled = false;
    this.siphonRefAngleBody = undefined; // will be recomputed on first step()
  }

  /**
   * Fire a siphon propulsion impulse.
   * 0° = jet exits along the vector from the aperture midpoint toward the
   * midpoint of CoM and CoB (hydrostatic stability axis direction in body frame).
   * Positive angles rotate CCW in body frame; 180° reverses the direction.
   * The animal accelerates in the *opposite* direction of the jet (Newton III).
   */
  applySiphonJet() {
    const m = this.model;
    const apex = m.centerline[m.centerline.length - 1];

    // Ensure siphonRefAngleBody is current; recompute if not yet set
    if (this.siphonRefAngleBody === undefined) {
      this._updateSiphonRefAngle(this.currentHydroState);
    }
    const baseAngle = (this.siphonRefAngleBody !== undefined)
      ? this.siphonRefAngleBody
      : apex.phi;

    // Siphon nozzle exits in direction (baseAngle + siphonAngleRad) in body space
    // The jet *exits* in that direction; reactive thrust is opposite (Newton III)
    const siphAngleRad = SIM_PARAMS.siphonAngle * Math.PI / 180;
    const jetLocalX = Math.cos(baseAngle + siphAngleRad);
    const jetLocalY = Math.sin(baseAngle + siphAngleRad);

    // Transform local body direction to world space
    const cosT = Math.cos(this.angle);
    const sinT = Math.sin(this.angle);
    const jetWorldX = jetLocalX * cosT - jetLocalY * sinT;
    const jetWorldY = jetLocalX * sinT + jetLocalY * cosT;

    // Reactive force opposes jet direction (Newton III); expressed in N, applied for ~0.15 s
    const F = SIM_PARAMS.siphonForce;
    this.siphonImpulse = {
      x: -jetWorldX * F,   // N, world space
      y: -jetWorldY * F,
      mag: F,
      // World position of aperture for torque calculation
      apX: apex.x - m.com.x, // body-local offset from CoM (mm)
      apY: apex.y - m.com.y,
      ttl: 0.15             // seconds of sustained thrust
    };
  }

  setDragStart(bodyPt, worldMouse) {
    this.isDragging = true;
    this.dragPickBody = { x: bodyPt.x, y: bodyPt.y };
    this.mouseWorld = { x: worldMouse.x, y: worldMouse.y };
  }

  setDragMove(worldMouse) {
    this.mouseWorld = { x: worldMouse.x, y: worldMouse.y };
  }

  setDragEnd() {
    this.isDragging = false;
    this.appliedForce = { x: 0, y: 0, mag: 0 };
    this.appliedTorque = 0;
  }

  step(dt) {
    const massKg = this.model.totalMass * 1e-3;
    const moiKgM2 = Math.max(1e-6, this.model.Izz_SI);

    // Evaluate hydrostatics at current pose
    const hydroState = this.hydro.evaluatePose(this.posX, this.posY, this.angle);

    // External user drag spring force
    let fUserX = 0;
    let fUserY = 0;
    let tauUser = 0;

    if (this.isDragging) {
      const cosT = Math.cos(this.angle);
      const sinT = Math.sin(this.angle);

      // Current world position of dragged pick point
      const pickWx = this.posX + (this.dragPickBody.x * cosT - this.dragPickBody.y * sinT);
      const pickWy = this.posY + (this.dragPickBody.x * sinT + this.dragPickBody.y * cosT);

      // Spring stiffness (N/mm -> N/m)
      const kSpring = 0.08; // N/mm
      const dx = this.mouseWorld.x - pickWx;
      const dy = this.mouseWorld.y - pickWy;

      fUserX = kSpring * dx;
      fUserY = kSpring * dy;

      const fMag = Math.hypot(fUserX, fUserY);
      this.appliedForce = { x: fUserX, y: fUserY, mag: fMag };

      // Torque from user force: r_arm x F_user (in N * m)
      const armX_m = (pickWx - this.posX) * 1e-3;
      const armY_m = (pickWy - this.posY) * 1e-3;
      tauUser = armX_m * fUserY - armY_m * fUserX;
      this.appliedTorque = tauUser;
    } else {
      this.appliedForce = { x: 0, y: 0, mag: 0 };
      this.appliedTorque = 0;
    }

    // Siphon jet impulse (spacebar-triggered, time-limited thrust)
    let fSiphX = 0;
    let fSiphY = 0;
    let tauSiph = 0;

    if (this.siphonImpulse && this.siphonImpulse.ttl > 0) {
      fSiphX = this.siphonImpulse.x;
      fSiphY = this.siphonImpulse.y;

      // Torque from siphon force applied at aperture (body-space offset → world)
      const cosT2 = Math.cos(this.angle);
      const sinT2 = Math.sin(this.angle);
      const armWx = (this.siphonImpulse.apX * cosT2 - this.siphonImpulse.apY * sinT2) * 1e-3; // m
      const armWy = (this.siphonImpulse.apX * sinT2 + this.siphonImpulse.apY * cosT2) * 1e-3;
      tauSiph = armWx * fSiphY - armWy * fSiphX; // N·m
      this.siphonImpulse.ttl -= dt;
    }

    // Hydrodynamic damping forces (water resistance)
    const subRatio = Math.min(1.0, hydroState.subVol_cm3 / this.model.totalOuterVolume);
    const dampingFactor = SIM_PARAMS.waterDamping ? 1.0 : 0.05;

    const cLin = 3.5 * massKg * dampingFactor * subRatio; // N * s / m
    const cAng = 4.0 * moiKgM2 * dampingFactor * subRatio; // N * m * s / rad

    // Net Forces (N)
    const fNetX = fUserX + fSiphX - cLin * (this.vx * 1e-3);
    const fNetY = (hydroState.buoyantForce - hydroState.gravityForce) + fUserY + fSiphY - cLin * (this.vy * 1e-3);
    const tauNet = hydroState.torque + tauUser + tauSiph - cAng * this.omega;

    // Linear and angular accelerations
    const ax_m = fNetX / massKg;
    const ay_m = fNetY / massKg;
    const alpha = tauNet / moiKgM2;

    // Numerical integration (semi-implicit Euler)
    this.vx += (ax_m * 1000) * dt; // mm/s
    this.vy += (ay_m * 1000) * dt;
    this.omega += alpha * dt;      // rad/s

    this.posX += this.vx * dt;
    this.posY += this.vy * dt;
    this.angle += this.omega * dt;

    // Restrict lateral drift smoothly
    this.posX *= 0.992;

    this.currentHydroState = hydroState;

    // Compute siphon reference direction in body-local frame:
    // 0° is defined as the vector FROM the aperture midpoint TOWARD the
    // midpoint between CoM (posX,posY) and CoB (currentHydroState.cob).
    // This aligns the neutral siphon axis with the hydrostatic stability couple.
    this._updateSiphonRefAngle(hydroState);
  }

  /**
   * Recompute siphonRefAngleBody from current pose + hydrostatic state.
   * Called after every physics step and before applySiphonJet().
   *
   * Result stored as this.siphonRefAngleBody (radians, body-local frame).
   * When the CoM == CoB (degenerate), falls back to aperture tangent direction.
   */
  _updateSiphonRefAngle(hydroState) {
    const m = this.model;
    const aperturePt = m.centerline[m.centerline.length - 1];
    const cosT = Math.cos(this.angle);
    const sinT = Math.sin(this.angle);

    // Aperture position in world space
    const apBodyX = aperturePt.x - m.com.x;
    const apBodyY = aperturePt.y - m.com.y;
    const apWorldX = this.posX + (apBodyX * cosT - apBodyY * sinT);
    const apWorldY = this.posY + (apBodyX * sinT + apBodyY * cosT);

    // CoM world position is (this.posX, this.posY)
    // CoB world position from hydroState
    const cobX = hydroState && hydroState.cob ? hydroState.cob.x : this.posX;
    const cobY = hydroState && hydroState.cob ? hydroState.cob.y : this.posY;

    // Midpoint between CoM and CoB in world space
    const midX = (this.posX + cobX) * 0.5;
    const midY = (this.posY + cobY) * 0.5;

    // Vector from aperture → midpoint, in world space
    const dxW = midX - apWorldX;
    const dyW = midY - apWorldY;
    const dMag = Math.hypot(dxW, dyW);

    if (dMag < 0.1) {
      // Degenerate: CoM == CoB (or aperture is at midpoint); fall back to aperture tangent
      this.siphonRefAngleBody = aperturePt.phi;
    } else {
      // Rotate world vector back into body frame (inverse rotation by -angle)
      const dxBody = dxW * cosT + dyW * sinT;
      const dyBody = -dxW * sinT + dyW * cosT;
      this.siphonRefAngleBody = Math.atan2(dyBody, dxBody);
    }
  }
}

/**
 * -------------------------------------------------------------
 * 4. High-Precision 2D Sagittal Cross-Section Renderer
 * -------------------------------------------------------------
 */
class Renderer2D {
  constructor(canvas, model, physics) {
    this.canvas = canvas;
    this.ctx = canvas.getContext("2d");
    this.model = model;
    this.physics = physics;

    this.zoom = 1.0;
    this.panX = 0;
    this.panY = 0;

    this.resize();
    window.addEventListener("resize", () => this.resize());
  }

  resize() {
    const rect = this.canvas.parentElement.getBoundingClientRect();
    this.width = rect.width;
    this.height = rect.height;
    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = this.width * dpr;
    this.canvas.height = this.height * dpr;
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  /**
   * World coordinates (mm) to Screen coordinates (px)
   */
  worldToScreen(wx, wy) {
    const scale = (Math.min(this.width, this.height) / (this.model.L * 1.35)) * this.zoom;
    const cx = this.width * 0.5 + this.panX;
    // Water surface (Y = 0) positioned at 90% the height of the panel (10% from the top)
    const cy = this.height * 0.10 + this.panY;
    return {
      x: cx + wx * scale,
      y: cy - wy * scale,
      scale
    };
  }

  /**
   * Screen coordinates (px) to World coordinates (mm)
   */
  screenToWorld(sx, sy) {
    const scale = (Math.min(this.width, this.height) / (this.model.L * 1.35)) * this.zoom;
    const cx = this.width * 0.5 + this.panX;
    const cy = this.height * 0.10 + this.panY;
    return {
      x: (sx - cx) / scale,
      y: (cy - sy) / scale
    };
  }

  render() {
    this.ctx.clearRect(0, 0, this.width, this.height);

    // Draw Water Medium Background
    this.drawWater();

    // Render Cephalopod Shell
    this.drawCephalopod();

    // Draw Hydrostatic & Biomechanical Overlays
    this.drawOverlays();
  }

  drawWater() {
    const screenWaterY = this.worldToScreen(0, 0).y;

    // Submerged water volume gradient
    const waterGrad = this.ctx.createLinearGradient(0, screenWaterY, 0, this.height);
    waterGrad.addColorStop(0, "rgba(14, 116, 144, 0.25)");
    waterGrad.addColorStop(1, "rgba(8, 51, 68, 0.75)");

    this.ctx.fillStyle = waterGrad;
    this.ctx.fillRect(0, screenWaterY, this.width, this.height - screenWaterY);

    // Water surface line with wave glow
    this.ctx.save();
    this.ctx.strokeStyle = "rgba(56, 189, 248, 0.85)";
    this.ctx.lineWidth = 2;
    this.ctx.shadowColor = "rgba(56, 189, 248, 0.6)";
    this.ctx.shadowBlur = 8;
    this.ctx.beginPath();
    this.ctx.moveTo(0, screenWaterY);
    this.ctx.lineTo(this.width, screenWaterY);
    this.ctx.stroke();

    // Waterline label
    this.ctx.fillStyle = "rgba(186, 230, 253, 0.85)";
    this.ctx.font = "11px ui-monospace, monospace";
    this.ctx.fillText("WATER SURFACE (Y = 0)", 16, screenWaterY - 6);
    this.ctx.restore();
  }

  drawCephalopod() {
    const p = this.physics;
    const m = this.model;
    const cosT = Math.cos(p.angle);
    const sinT = Math.sin(p.angle);
    const com = m.com;

    this.ctx.save();

    // Transform points from local shell space to screen
    const transformPt = (lx, ly) => {
      const rx = lx - com.x;
      const ry = ly - com.y;
      const wx = p.posX + (rx * cosT - ry * sinT);
      const wy = p.posY + (rx * sinT + ry * cosT);
      return this.worldToScreen(wx, wy);
    };

    // 1. Phragmocone Chambers (Air Fill)
    this.ctx.beginPath();
    for (let i = 0; i < m.innerVentral.length; i++) {
      if (m.centerline[i].s > (m.L - m.L_bc)) break;
      const pt = transformPt(m.innerVentral[i].x, m.innerVentral[i].y);
      if (i === 0) this.ctx.moveTo(pt.x, pt.y);
      else this.ctx.lineTo(pt.x, pt.y);
    }
    for (let i = m.innerDorsal.length - 1; i >= 0; i--) {
      if (m.centerline[i].s > (m.L - m.L_bc)) continue;
      const pt = transformPt(m.innerDorsal[i].x, m.innerDorsal[i].y);
      this.ctx.lineTo(pt.x, pt.y);
    }
    this.ctx.closePath();
    this.ctx.fillStyle = "rgba(224, 242, 254, 0.08)"; // Gas-filled camerae
    this.ctx.fill();

    // 2. Living Chamber Soft Visceral Body (Supplementary Weight)
    this.ctx.beginPath();
    let started = false;
    for (let i = 0; i < m.innerVentral.length; i++) {
      if (m.centerline[i].s >= (m.L - m.L_bc)) {
        const pt = transformPt(m.innerVentral[i].x, m.innerVentral[i].y);
        if (!started) {
          this.ctx.moveTo(pt.x, pt.y);
          started = true;
        } else {
          this.ctx.lineTo(pt.x, pt.y);
        }
      }
    }
    for (let i = m.innerDorsal.length - 1; i >= 0; i--) {
      if (m.centerline[i].s >= (m.L - m.L_bc)) {
        const pt = transformPt(m.innerDorsal[i].x, m.innerDorsal[i].y);
        this.ctx.lineTo(pt.x, pt.y);
      }
    }
    this.ctx.closePath();
    const bodyGrad = this.ctx.createLinearGradient(
      this.width * 0.4, this.height * 0.4,
      this.width * 0.6, this.height * 0.6
    );
    bodyGrad.addColorStop(0, "rgba(225, 29, 72, 0.4)");
    bodyGrad.addColorStop(1, "rgba(190, 24, 93, 0.6)");
    this.ctx.fillStyle = bodyGrad;
    this.ctx.fill();

    // 3. Conch Mineralized Wall (Calcite/Aragonite)
    this.ctx.beginPath();
    // Outer boundary (Ventral -> Aperture -> Dorsal -> Apex)
    for (let i = 0; i < m.ventralWall.length; i++) {
      const pt = transformPt(m.ventralWall[i].x, m.ventralWall[i].y);
      if (i === 0) this.ctx.moveTo(pt.x, pt.y);
      else this.ctx.lineTo(pt.x, pt.y);
    }
    for (let i = m.dorsalWall.length - 1; i >= 0; i--) {
      const pt = transformPt(m.dorsalWall[i].x, m.dorsalWall[i].y);
      this.ctx.lineTo(pt.x, pt.y);
    }
    this.ctx.closePath();

    // Cutout inner boundary using even-odd fill rule
    for (let i = 0; i < m.innerVentral.length; i++) {
      const pt = transformPt(m.innerVentral[i].x, m.innerVentral[i].y);
      if (i === 0) this.ctx.moveTo(pt.x, pt.y);
      else this.ctx.lineTo(pt.x, pt.y);
    }
    for (let i = m.innerDorsal.length - 1; i >= 0; i--) {
      const pt = transformPt(m.innerDorsal[i].x, m.innerDorsal[i].y);
      this.ctx.lineTo(pt.x, pt.y);
    }
    this.ctx.closePath();

    this.ctx.fillStyle = "rgba(217, 180, 130, 0.88)"; // Aragonite / Calcite shell
    this.ctx.strokeStyle = "#854d0e";
    this.ctx.lineWidth = 1.2;
    this.ctx.fill("evenodd");
    this.ctx.stroke();

    // 4. Mineralized Septa (transverse dividing walls)
    this.ctx.strokeStyle = "rgba(217, 180, 130, 0.95)";
    this.ctx.lineWidth = Math.max(1.0, m.t_septum * 0.8);

    for (const sep of m.septa) {
      // Dome-shaped curved septal arc
      const vPt = transformPt(sep.x + sep.r_in * sep.nx, sep.y + sep.r_in * sep.ny);
      const dPt = transformPt(sep.x - sep.r_in * sep.nx, sep.y - sep.r_in * sep.ny);
      // Septal curvature convexity pointing toward apex
      const cPt = transformPt(sep.x - sep.r_in * 0.3 * Math.cos(sep.phi), sep.y - sep.r_in * 0.3 * Math.sin(sep.phi));

      this.ctx.beginPath();
      this.ctx.moveTo(vPt.x, vPt.y);
      this.ctx.quadraticCurveTo(cPt.x, cPt.y, dPt.x, dPt.y);
      this.ctx.stroke();
    }

    // 5. Siphuncular Conical Tube (running through phragmocone)
    if (m.siphuncle && m.siphuncle.length > 1) {
      this.ctx.beginPath();
      // Ventral surface of siphuncle cone from apex to phragmocone opening
      for (let i = 0; i < m.siphuncle.length; i++) {
        const pt = transformPt(m.siphuncle[i].vx, m.siphuncle[i].vy);
        if (i === 0) this.ctx.moveTo(pt.x, pt.y);
        else this.ctx.lineTo(pt.x, pt.y);
      }
      // Dorsal surface of siphuncle cone back to apex
      for (let i = m.siphuncle.length - 1; i >= 0; i--) {
        const pt = transformPt(m.siphuncle[i].dx, m.siphuncle[i].dy);
        this.ctx.lineTo(pt.x, pt.y);
      }
      this.ctx.closePath();
      this.ctx.fillStyle = "rgba(217, 119, 6, 0.55)"; // Mineralized siphuncular cone
      this.ctx.strokeStyle = "#9a3412";
      this.ctx.lineWidth = 1.4;
      this.ctx.fill();
      this.ctx.stroke();

      // Central siphuncular vascular lumen line
      this.ctx.beginPath();
      for (let i = 0; i < m.siphuncle.length; i++) {
        const pt = transformPt(m.siphuncle[i].cx, m.siphuncle[i].cy);
        if (i === 0) this.ctx.moveTo(pt.x, pt.y);
        else this.ctx.lineTo(pt.x, pt.y);
      }
      this.ctx.strokeStyle = "rgba(254, 215, 170, 0.85)";
      this.ctx.lineWidth = 1.0;
      this.ctx.stroke();
    }

    this.ctx.restore();
  }

  drawOverlays() {
    const p = this.physics;
    const m = this.model;
    const cosT = Math.cos(p.angle);
    const sinT = Math.sin(p.angle);

    // Screen position of Center of Mass (CoM / CG)
    const comScr = this.worldToScreen(p.posX, p.posY);

    // 1. Moment of Inertia (Principal Inertia Ellipse and Axes)
    if (SIM_PARAMS.showMoI) {
      this.ctx.save();
      this.ctx.translate(comScr.x, comScr.y);
      // Rotate by body angle + principal inertia angle
      this.ctx.rotate(-(p.angle + m.principalAngle));

      const scale = comScr.scale;
      const rad1 = m.k1 * scale; // Major radius of gyration
      const rad2 = m.k2 * scale; // Minor radius of gyration

      // Draw inertia ellipse
      this.ctx.beginPath();
      this.ctx.ellipse(0, 0, Math.max(2, rad1), Math.max(2, rad2), 0, 0, 2 * Math.PI);
      this.ctx.strokeStyle = "rgba(245, 158, 11, 0.75)";
      this.ctx.lineWidth = 1.8;
      this.ctx.setLineDash([5, 4]);
      this.ctx.fillStyle = "rgba(245, 158, 11, 0.06)";
      this.ctx.fill();
      this.ctx.stroke();

      // Principal axes
      this.ctx.beginPath();
      this.ctx.moveTo(-rad1 * 1.2, 0);
      this.ctx.lineTo(rad1 * 1.2, 0);
      this.ctx.moveTo(0, -rad2 * 1.2);
      this.ctx.lineTo(0, rad2 * 1.2);
      this.ctx.strokeStyle = "rgba(245, 158, 11, 0.4)";
      this.ctx.stroke();

      this.ctx.restore();
    }

    // 2. Center of Buoyancy (CoB / B)
    if (SIM_PARAMS.showCoB && p.currentHydroState) {
      const cobScr = this.worldToScreen(p.currentHydroState.cob.x, p.currentHydroState.cob.y);

      // Metacentric separation line (B - G couple)
      this.ctx.save();
      this.ctx.beginPath();
      this.ctx.moveTo(comScr.x, comScr.y);
      this.ctx.lineTo(cobScr.x, cobScr.y);
      this.ctx.strokeStyle = "rgba(255, 255, 255, 0.5)";
      this.ctx.lineWidth = 1.5;
      this.ctx.setLineDash([3, 3]);
      this.ctx.stroke();

      // Buoyancy force vector arrow upward
      const fBuoy = p.currentHydroState.buoyantForce;
      const arrowLen = Math.min(80, fBuoy * 3.5);
      this.drawArrow(cobScr.x, cobScr.y, cobScr.x, cobScr.y - arrowLen, "#00d2ff", "F_B");

      // CoB Target Marker
      this.drawTargetMarker(cobScr.x, cobScr.y, "#00d2ff", "CoB (B)");
      this.ctx.restore();
    }

    // 3. Center of Mass (CoM / G)
    if (SIM_PARAMS.showCoM) {
      this.ctx.save();
      // Gravity force arrow downward
      const fGrav = (m.totalMass * 1e-3) * GRAVITY;
      const arrowLen = Math.min(80, fGrav * 3.5);
      this.drawArrow(comScr.x, comScr.y, comScr.x, comScr.y + arrowLen, "#ef4444", "F_g");

      // CoM Target Marker
      this.drawTargetMarker(comScr.x, comScr.y, "#ef4444", "CoM (G)");
      this.ctx.restore();
    }

    // 4. Interactive User Drag Force Vector
    if (p.isDragging && SIM_PARAMS.showForces) {
      const pickWx = p.posX + (p.dragPickBody.x * cosT - p.dragPickBody.y * sinT);
      const pickWy = p.posY + (p.dragPickBody.x * sinT + p.dragPickBody.y * cosT);
      const pickScr = this.worldToScreen(pickWx, pickWy);
      const mouseScr = this.worldToScreen(p.mouseWorld.x, p.mouseWorld.y);

      this.ctx.save();
      // Spring line
      this.ctx.beginPath();
      this.ctx.moveTo(pickScr.x, pickScr.y);
      this.ctx.lineTo(mouseScr.x, mouseScr.y);
      this.ctx.strokeStyle = "#eab308";
      this.ctx.lineWidth = 2.5;
      this.ctx.shadowColor = "rgba(234, 179, 8, 0.6)";
      this.ctx.shadowBlur = 6;
      this.ctx.stroke();

      // Pick point handle
      this.ctx.beginPath();
      this.ctx.arc(pickScr.x, pickScr.y, 6, 0, 2 * Math.PI);
      this.ctx.fillStyle = "#eab308";
      this.ctx.fill();

      // Force Vector Arrow & Label
      const fN = p.appliedForce.mag;
      const fGf = fN * 101.97;
      this.drawArrow(pickScr.x, pickScr.y, mouseScr.x, mouseScr.y, "#ffeb3b", "");

      this.ctx.fillStyle = "#ffeb3b";
      this.ctx.font = "bold 12px ui-monospace, monospace";
      this.ctx.fillText(`F = ${fN.toFixed(2)} N (${fGf.toFixed(0)} gf)`, mouseScr.x + 10, mouseScr.y - 10);

      this.ctx.restore();
    }

    // 5. Hyponome (Siphon) — nozzle arrow at aperture + firing flash
    {
      const aperturePt = m.centerline[m.centerline.length - 1]; // last centerline point = aperture

      // Base angle: vector from aperture → CoM-CoB midpoint in body frame (or aperture tangent if undefined)
      const baseAngle = (p.siphonRefAngleBody !== undefined) ? p.siphonRefAngleBody : aperturePt.phi;

      // Aperture local pos → world
      const apWorldX = p.posX + ((aperturePt.x - m.com.x) * cosT - (aperturePt.y - m.com.y) * sinT);
      const apWorldY = p.posY + ((aperturePt.x - m.com.x) * sinT + (aperturePt.y - m.com.y) * cosT);
      const apScr = this.worldToScreen(apWorldX, apWorldY);

      // Jet exit direction in world space (the direction water squirts OUT)
      const siphAngleRad = SIM_PARAMS.siphonAngle * Math.PI / 180;
      const jetLocalX = Math.cos(baseAngle + siphAngleRad);
      const jetLocalY = Math.sin(baseAngle + siphAngleRad);
      const jetWorldX = jetLocalX * cosT - jetLocalY * sinT;
      const jetWorldY = jetLocalX * sinT + jetLocalY * cosT;

      // Nozzle arrow length proportional to force (min 40 px, max 100 px)
      const arrowPx = Math.min(100, Math.max(40, SIM_PARAMS.siphonForce * 35));

      // Screen endpoint: jet exit direction (arrow points the way water goes)
      const nozzleTipX = apScr.x + jetWorldX * arrowPx;
      const nozzleTipY = apScr.y - jetWorldY * arrowPx;   // y is inverted in screen space

      // Animal reaction direction (opposite of jet)
      const reactTipX = apScr.x - jetWorldX * arrowPx;
      const reactTipY = apScr.y + jetWorldY * arrowPx;

      this.ctx.save();

      // 0° reference tick — faint dotted line toward CoM-CoB midpoint (the neutral axis)
      const refLocalX = Math.cos(baseAngle);
      const refLocalY = Math.sin(baseAngle);
      const refWorldX = refLocalX * cosT - refLocalY * sinT;
      const refWorldY = refLocalX * sinT + refLocalY * cosT;
      const refTickLen = arrowPx * 0.55;
      this.ctx.beginPath();
      this.ctx.setLineDash([2, 5]);
      this.ctx.moveTo(apScr.x, apScr.y);
      this.ctx.lineTo(apScr.x + refWorldX * refTickLen, apScr.y - refWorldY * refTickLen);
      this.ctx.strokeStyle = "rgba(167, 243, 208, 0.30)";
      this.ctx.lineWidth = 1.2;
      this.ctx.stroke();
      this.ctx.setLineDash([]);

      // Dashed nozzle guide line (jet exit)
      this.ctx.beginPath();
      this.ctx.setLineDash([6, 4]);
      this.ctx.moveTo(apScr.x, apScr.y);
      this.ctx.lineTo(nozzleTipX, nozzleTipY);
      this.ctx.strokeStyle = "rgba(167, 243, 208, 0.65)";
      this.ctx.lineWidth = 1.8;
      this.ctx.stroke();
      this.ctx.setLineDash([]);

      // Siphon nozzle funnel at aperture (small filled circle)
      this.ctx.beginPath();
      this.ctx.arc(apScr.x, apScr.y, 5, 0, 2 * Math.PI);
      this.ctx.fillStyle = "#34d399";
      this.ctx.fill();

      // Jet exit label
      this.ctx.fillStyle = "#34d399";
      this.ctx.font = "11px ui-monospace, monospace";
      this.ctx.fillText(`Jet (${SIM_PARAMS.siphonAngle.toFixed(0)}°)`, nozzleTipX + 6, nozzleTipY - 4);

      // Reaction thrust arrow (direction animal will move) — solid green
      this.drawArrow(apScr.x, apScr.y, reactTipX, reactTipY, "#10b981", `F_jet = ${SIM_PARAMS.siphonForce.toFixed(1)} N`);

      // If siphon jet is currently firing, draw burst flash
      if (p.siphonImpulse && p.siphonImpulse.ttl > 0) {
        const flashAlpha = Math.min(1.0, p.siphonImpulse.ttl / 0.15);
        const burstLen = arrowPx * 1.6;
        const burstTipX = apScr.x + jetWorldX * burstLen;
        const burstTipY = apScr.y - jetWorldY * burstLen;

        // Gradient jet plume
        const grad = this.ctx.createLinearGradient(apScr.x, apScr.y, burstTipX, burstTipY);
        grad.addColorStop(0, `rgba(52, 211, 153, ${flashAlpha * 0.95})`);
        grad.addColorStop(1, `rgba(52, 211, 153, 0)`);

        this.ctx.beginPath();
        this.ctx.moveTo(apScr.x, apScr.y);
        this.ctx.lineTo(burstTipX, burstTipY);
        this.ctx.strokeStyle = grad;
        this.ctx.lineWidth = 6 * flashAlpha;
        this.ctx.shadowColor = "rgba(52, 211, 153, 0.8)";
        this.ctx.shadowBlur = 14 * flashAlpha;
        this.ctx.stroke();
        this.ctx.shadowBlur = 0;
      }

      this.ctx.restore();
    }
  }

  drawTargetMarker(x, y, color, label) {
    this.ctx.save();
    this.ctx.strokeStyle = color;
    this.ctx.fillStyle = color;
    this.ctx.lineWidth = 2;

    // Outer circle
    this.ctx.beginPath();
    this.ctx.arc(x, y, 7, 0, 2 * Math.PI);
    this.ctx.stroke();

    // Inner dot
    this.ctx.beginPath();
    this.ctx.arc(x, y, 2.5, 0, 2 * Math.PI);
    this.ctx.fill();

    // Crosshairs
    this.ctx.beginPath();
    this.ctx.moveTo(x - 10, y);
    this.ctx.lineTo(x + 10, y);
    this.ctx.moveTo(x, y - 10);
    this.ctx.lineTo(x, y + 10);
    this.ctx.stroke();

    // Label
    this.ctx.font = "bold 11px ui-monospace, monospace";
    this.ctx.fillText(label, x + 12, y + 4);
    this.ctx.restore();
  }

  drawArrow(fromX, fromY, toX, toY, color, label) {
    const headLen = 9;
    const dx = toX - fromX;
    const dy = toY - fromY;
    const angle = Math.atan2(dy, dx);

    this.ctx.save();
    this.ctx.strokeStyle = color;
    this.ctx.fillStyle = color;
    this.ctx.lineWidth = 2;

    this.ctx.beginPath();
    this.ctx.moveTo(fromX, fromY);
    this.ctx.lineTo(toX, toY);
    this.ctx.stroke();

    this.ctx.beginPath();
    this.ctx.moveTo(toX, toY);
    this.ctx.lineTo(toX - headLen * Math.cos(angle - Math.PI / 6), toY - headLen * Math.sin(angle - Math.PI / 6));
    this.ctx.lineTo(toX - headLen * Math.cos(angle + Math.PI / 6), toY - headLen * Math.sin(angle + Math.PI / 6));
    this.ctx.closePath();
    this.ctx.fill();

    if (label) {
      this.ctx.font = "bold 11px ui-monospace, monospace";
      this.ctx.fillText(label, toX + 6, toY + 4);
    }
    this.ctx.restore();
  }

  /**
   * Hit test to find if a world point is within the cephalopod shell
   */
  hitTest(worldX, worldY) {
    const p = this.physics;
    const m = this.model;
    const com = m.com;

    // Transform world point to shell local coordinates
    const cosT = Math.cos(-p.angle);
    const sinT = Math.sin(-p.angle);
    const dx = worldX - p.posX;
    const dy = worldY - p.posY;
    const lx = (dx * cosT - dy * sinT) + com.x;
    const ly = (dx * sinT + dy * cosT) + com.y;

    // Find nearest centerline point
    let minD2 = Infinity;
    let nearestSlice = null;

    for (const slice of m.slices) {
      const d2 = (slice.x - lx) * (slice.x - lx) + (slice.y - ly) * (slice.y - ly);
      if (d2 < minD2) {
        minD2 = d2;
        nearestSlice = slice;
      }
    }

    if (nearestSlice && minD2 <= (nearestSlice.r_out * 1.3) * (nearestSlice.r_out * 1.3)) {
      return {
        hit: true,
        bodyPt: { x: lx - com.x, y: ly - com.y }
      };
    }

    return { hit: false };
  }
}

/**
 * -------------------------------------------------------------
 * 5. Three.js WebGL 3D Perspective Visualizer
 * -------------------------------------------------------------
 */
class Renderer3D {
  constructor(container, model, physics) {
    this.container = container;
    this.model = model;
    this.physics = physics;
    this.initialized = false;
  }

  init() {
    if (this.initialized) return;

    const width = this.container.clientWidth || 800;
    const height = this.container.clientHeight || 600;

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x080d14);

    this.camera = new THREE.PerspectiveCamera(45, width / height, 1, 5000);
    this.camera.position.set(0, 150, 750);

    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setSize(width, height);
    this.renderer.setPixelRatio(window.devicePixelRatio || 1);
    this.renderer.shadowMap.enabled = true;
    this.container.appendChild(this.renderer.domElement);

    if (typeof THREE.OrbitControls !== "undefined") {
      this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
      this.controls.enableDamping = true;
      this.controls.dampingFactor = 0.05;
    }

    // Lighting
    const ambLight = new THREE.AmbientLight(0xffffff, 0.45);
    this.scene.add(ambLight);

    const dirLight = new THREE.DirectionalLight(0xffffff, 0.85);
    dirLight.position.set(300, 600, 400);
    this.scene.add(dirLight);

    // Water plane
    const waterGeom = new THREE.PlaneGeometry(3000, 3000, 32, 32);
    const waterMat = new THREE.MeshStandardMaterial({
      color: 0x0284c7,
      transparent: true,
      opacity: 0.35,
      roughness: 0.1,
      metalness: 0.2,
      side: THREE.DoubleSide
    });
    this.waterMesh = new THREE.Mesh(waterGeom, waterMat);
    this.waterMesh.rotation.x = -Math.PI / 2;
    this.scene.add(this.waterMesh);

    // Cephalopod 3D Root Group
    this.shellGroup = new THREE.Group();
    this.scene.add(this.shellGroup);

    this.buildShellMesh();
    this.initialized = true;

    window.addEventListener("resize", () => {
      const w = this.container.clientWidth;
      const h = this.container.clientHeight;
      this.camera.aspect = w / h;
      this.camera.updateProjectionMatrix();
      this.renderer.setSize(w, h);
    });
  }

  buildShellMesh() {
    // Remove old meshes
    while (this.shellGroup.children.length > 0) {
      this.shellGroup.remove(this.shellGroup.children[0]);
    }

    const m = this.model;
    const com = m.com;

    // Create 3D sweep geometry along centerline
    const curvePoints = [];
    const radii = [];

    for (let i = 0; i <= m.numSlices; i += 2) {
      const pt = m.centerline[i];
      // Centered relative to CoM
      curvePoints.push(new THREE.Vector3(pt.x - com.x, pt.y - com.y, 0));
      radii.push(pt.r_out);
    }

    const curve = new THREE.CatmullRomCurve3(curvePoints);
    const tubeGeom = new THREE.TubeGeometry(curve, 100, (m.D_ap / 4), 18, false);

    // Shell material: pearly nacre / calcite with ribbed normal texture
    const shellMat = new THREE.MeshStandardMaterial({
      color: 0xd4b886,
      roughness: 0.45,
      metalness: 0.15,
      side: THREE.DoubleSide
    });

    const shellMesh = new THREE.Mesh(tubeGeom, shellMat);
    this.shellGroup.add(shellMesh);

    // Markers for CoM and CoB
    const comGeom = new THREE.SphereGeometry(6, 16, 16);
    const comMat = new THREE.MeshBasicMaterial({ color: 0xef4444 });
    const comMesh = new THREE.Mesh(comGeom, comMat);
    this.shellGroup.add(comMesh);

    // Inertia Ellipsoid Wireframe
    const ellipGeom = new THREE.SphereGeometry(1, 16, 12);
    const ellipMat = new THREE.MeshBasicMaterial({
      color: 0xf59e0b,
      wireframe: true,
      transparent: true,
      opacity: 0.4
    });
    this.ellipMesh = new THREE.Mesh(ellipGeom, ellipMat);
    this.ellipMesh.scale.set(m.k1, m.k2, (m.k1 + m.k2) * 0.5);
    this.shellGroup.add(this.ellipMesh);

    // 3D Siphuncular Cone
    if (m.siphuncle && m.siphuncle.length > 2) {
      const siphCurvePoints = [];
      for (let i = 0; i < m.siphuncle.length; i += 2) {
        const pt = m.siphuncle[i];
        siphCurvePoints.push(new THREE.Vector3(pt.cx - com.x, pt.cy - com.y, 0));
      }
      if (siphCurvePoints.length > 1) {
        const siphCurve = new THREE.CatmullRomCurve3(siphCurvePoints);
        const siphTubeGeom = new THREE.TubeGeometry(siphCurve, 50, Math.max(0.6, m.siphDiam * 0.35), 10, false);
        const siphMat = new THREE.MeshStandardMaterial({
          color: 0xb45309,
          roughness: 0.35,
          metalness: 0.2
        });
        const siphMesh = new THREE.Mesh(siphTubeGeom, siphMat);
        this.shellGroup.add(siphMesh);
      }
    }
  }

  render() {
    if (!this.initialized) return;

    if (this.controls) this.controls.update();

    // Sync shell pose with physics simulation
    const p = this.physics;
    this.shellGroup.position.set(p.posX, p.posY, 0);
    this.shellGroup.rotation.z = p.angle;

    this.renderer.render(this.scene, this.camera);
  }
}

/**
 * -------------------------------------------------------------
 * 6. Main UI Controller & Simulation Orchestration
 * -------------------------------------------------------------
 */
class SimulationApp {
  constructor() {
    this.model = new CephalopodModel(SIM_PARAMS);
    this.hydro = new HydrostaticSolver(this.model);
    this.physics = new PhysicsEngine(this.model, this.hydro);

    this.canvas2D = document.getElementById("sim-canvas-2d");
    this.container3D = document.getElementById("sim-container-3d");
    this.renderer2D = new Renderer2D(this.canvas2D, this.model, this.physics);
    this.renderer3D = new Renderer3D(this.container3D, this.model, this.physics);

    this.is3DMode = false;
    this.isPaused = false;

    this.bindControls();
    this.bindMouseInteractions();
    this.updateHUD();

    this.lastTime = performance.now();
    this.animate = this.animate.bind(this);
    requestAnimationFrame(this.animate);
  }

  bindControls() {
    // Presets
    const presetSelect = document.getElementById("preset-select");
    presetSelect.addEventListener("change", (e) => {
      const key = e.target.value;
      if (CEPHALOPOD_PRESETS[key]) {
        const p = CEPHALOPOD_PRESETS[key];
        document.getElementById("input-length").value = p.shellLength;
        document.getElementById("input-curvature").value = p.curvature;
        document.getElementById("input-aperture-diam").value = p.apertureDiam;
        document.getElementById("input-living-length").value = p.livingChamberLength;
        document.getElementById("input-shell-thickness").value = p.shellThickness;
        document.getElementById("input-septa-thickness").value = p.septaThickness;
        document.getElementById("input-chamber-dist").value = p.chamberDist;
        if (p.siphuncleDiam !== undefined) document.getElementById("input-siphuncle-diam").value = p.siphuncleDiam;
        if (p.siphunclePos !== undefined) document.getElementById("input-siphuncle-pos").value = p.siphunclePos;
        if (p.siphuncleDensity !== undefined) document.getElementById("input-siphuncle-density").value = p.siphuncleDensity;
        document.getElementById("input-shell-density").value = p.shellDensity;
        document.getElementById("input-chamber-density").value = p.chamberGasDensity !== undefined ? p.chamberGasDensity : 0.0012;
        document.getElementById("input-weight-factor").value = p.weightFactor;
        document.getElementById("input-cameral-fluid").value = p.cameralFluid;
        document.getElementById("select-water-density").value = p.waterDensity;
        this.onParamsChanged();
      }
    });

    // Parameter sliders
    const sliderBindings = [
      { id: "input-length", param: "shellLength", disp: "val-length", unit: " mm" },
      { id: "input-curvature", param: "curvature", disp: "val-curvature", unit: "", fmt: (v) => {
        if (v <= 0.05) return `${v.toFixed(2)} (Orthocone)`;
        if (v <= 0.35) return `${v.toFixed(2)} (Cyrtocone)`;
        if (v <= 0.65) return `${v.toFixed(2)} (Gyrocone)`;
        return `${v.toFixed(2)} (Nautilocone)`;
      }},
      { id: "input-aperture-diam", param: "apertureDiam", disp: "val-aperture-diam", unit: " mm" },
      { id: "input-living-length", param: "livingChamberLength", disp: "val-living-length", unit: " mm" },
      { id: "input-shell-thickness", param: "shellThickness", disp: "val-shell-thickness", unit: " mm" },
      { id: "input-septa-thickness", param: "septaThickness", disp: "val-septa-thickness", unit: " mm" },
      { id: "input-chamber-dist", param: "chamberDist", disp: "val-chamber-dist", unit: " mm" },
      { id: "input-siphuncle-diam", param: "siphuncleDiam", disp: "val-siphuncle-diam", unit: " mm" },
      { id: "input-siphuncle-pos", param: "siphunclePos", disp: "val-siphuncle-pos", unit: "", fmt: (v) => {
        if (v <= 5) return `${v}% (Ventral)`;
        if (v >= 45 && v <= 55) return `${v}% (Center)`;
        if (v >= 95) return `${v}% (Dorsal)`;
        return `${v}%`;
      }},
      { id: "input-siphuncle-density", param: "siphuncleDensity", disp: "val-siphuncle-density", unit: " g/cm³" },
      { id: "input-shell-density", param: "shellDensity", disp: "val-shell-density", unit: " g/cm³" },
      { id: "input-chamber-density", param: "chamberGasDensity", disp: "val-chamber-density", unit: "", fmt: (v) => {
        if (v <= 0.0005) return `${v.toFixed(4)} g/cm³ (Near Vacuum)`;
        if (v <= 0.0015) return `${v.toFixed(4)} g/cm³ (Standard Air)`;
        if (v <= 0.003) return `${v.toFixed(4)} g/cm³ (CO₂ Gas)`;
        if (v >= 0.95 && v <= 1.05) return `${v.toFixed(3)} g/cm³ (Water/Brine)`;
        return `${v.toFixed(4)} g/cm³`;
      }},
      { id: "input-weight-factor", param: "weightFactor", disp: "val-weight-factor", unit: "", fmt: (v) => `${v.toFixed(2)} (${(1.05 * v).toFixed(2)} g/cm³)` },
      { id: "input-cameral-fluid", param: "cameralFluid", disp: "val-cameral-fluid", unit: " %" }
    ];

    sliderBindings.forEach((b) => {
      const el = document.getElementById(b.id);
      const disp = document.getElementById(b.disp);
      el.addEventListener("input", (e) => {
        const val = parseFloat(e.target.value);
        SIM_PARAMS[b.param] = val;
        disp.textContent = b.fmt ? b.fmt(val) : `${val}${b.unit}`;
        presetSelect.value = "custom";
        this.onParamsChanged();
      });
    });

    // Siphuncle density synchronization with shell density
    const chkSyncSiph = document.getElementById("chk-sync-siph-density");
    const inputShellDensity = document.getElementById("input-shell-density");
    const inputSiphDensity = document.getElementById("input-siphuncle-density");
    const dispSiphDensity = document.getElementById("val-siphuncle-density");

    inputShellDensity.addEventListener("input", () => {
      if (chkSyncSiph && chkSyncSiph.checked) {
        inputSiphDensity.value = SIM_PARAMS.shellDensity;
        SIM_PARAMS.siphuncleDensity = SIM_PARAMS.shellDensity;
        dispSiphDensity.textContent = `${SIM_PARAMS.shellDensity.toFixed(2)} g/cm³`;
      }
    });

    inputSiphDensity.addEventListener("input", () => {
      if (chkSyncSiph) chkSyncSiph.checked = false;
    });

    if (chkSyncSiph) {
      chkSyncSiph.addEventListener("change", (e) => {
        if (e.target.checked) {
          inputSiphDensity.value = SIM_PARAMS.shellDensity;
          SIM_PARAMS.siphuncleDensity = SIM_PARAMS.shellDensity;
          dispSiphDensity.textContent = `${SIM_PARAMS.shellDensity.toFixed(2)} g/cm³`;
          this.onParamsChanged();
        }
      });
    }

    document.getElementById("select-water-density").addEventListener("change", (e) => {
      SIM_PARAMS.waterDensity = parseFloat(e.target.value);
      this.onParamsChanged();
    });

    // Checkboxes
    document.getElementById("chk-show-com").addEventListener("change", (e) => SIM_PARAMS.showCoM = e.target.checked);
    document.getElementById("chk-show-cob").addEventListener("change", (e) => SIM_PARAMS.showCoB = e.target.checked);
    document.getElementById("chk-show-moi").addEventListener("change", (e) => SIM_PARAMS.showMoI = e.target.checked);
    document.getElementById("chk-show-forces").addEventListener("change", (e) => SIM_PARAMS.showForces = e.target.checked);
    document.getElementById("chk-water-damping").addEventListener("change", (e) => SIM_PARAMS.waterDamping = e.target.checked);

    // Siphon / Hyponome controls (angle slider + force input)
    const siphonAngleEl = document.getElementById("input-siphon-angle");
    const siphonAngleDisp = document.getElementById("val-siphon-angle");
    if (siphonAngleEl) {
      siphonAngleEl.addEventListener("input", (e) => {
        SIM_PARAMS.siphonAngle = parseFloat(e.target.value);
        const a = SIM_PARAMS.siphonAngle;
        let label = `${a.toFixed(0)}°`;
        if (a <= 5) label += " (Forward)";
        else if (a >= 85 && a <= 95) label += " (Perpendicular)";
        else if (a >= 175) label += " (Backward)";
        siphonAngleDisp.textContent = label;
      });
    }

    const siphonForceEl = document.getElementById("input-siphon-force");
    const siphonForceDisp = document.getElementById("val-siphon-force");
    if (siphonForceEl) {
      siphonForceEl.addEventListener("input", (e) => {
        SIM_PARAMS.siphonForce = parseFloat(e.target.value);
        siphonForceDisp.textContent = `${SIM_PARAMS.siphonForce.toFixed(1)} N`;
      });
    }

    // Spacebar → fire siphon jet impulse
    window.addEventListener("keydown", (e) => {
      if (e.code === "Space" && !e.repeat) {
        e.preventDefault(); // prevent page scroll
        if (!this.isPaused) {
          this.physics.applySiphonJet();
        }
      }
    });

    // Buttons
    // 1. Pause / Resume button
    const btnPause = document.getElementById("btn-pause-phys");
    btnPause.addEventListener("click", () => {
      this.isPaused = !this.isPaused;
      btnPause.textContent = this.isPaused ? "Resume" : "Pause";
      btnPause.classList.toggle("btn-primary", this.isPaused);
      btnPause.classList.toggle("btn-secondary", !this.isPaused);
    });

    // 2. Reset All to Defaults button
    document.getElementById("btn-reset-all").addEventListener("click", () => {
      const defaultPreset = CEPHALOPOD_PRESETS.orthoceras;
      presetSelect.value = "orthoceras";
      document.getElementById("input-length").value = defaultPreset.shellLength;
      document.getElementById("input-curvature").value = defaultPreset.curvature;
      document.getElementById("input-aperture-diam").value = defaultPreset.apertureDiam;
      document.getElementById("input-living-length").value = defaultPreset.livingChamberLength;
      document.getElementById("input-shell-thickness").value = defaultPreset.shellThickness;
      document.getElementById("input-septa-thickness").value = defaultPreset.septaThickness;
      document.getElementById("input-chamber-dist").value = defaultPreset.chamberDist;
      document.getElementById("input-siphuncle-diam").value = defaultPreset.siphuncleDiam;
      document.getElementById("input-siphuncle-pos").value = defaultPreset.siphunclePos;
      document.getElementById("input-siphuncle-density").value = defaultPreset.siphuncleDensity;
      if (chkSyncSiph) chkSyncSiph.checked = true;
      document.getElementById("input-shell-density").value = defaultPreset.shellDensity;
      document.getElementById("input-chamber-density").value = defaultPreset.chamberGasDensity;
      document.getElementById("input-weight-factor").value = defaultPreset.weightFactor;
      document.getElementById("input-cameral-fluid").value = defaultPreset.cameralFluid;
      document.getElementById("select-water-density").value = defaultPreset.waterDensity;

      // Update slider displays
      sliderBindings.forEach((b) => {
        const el = document.getElementById(b.id);
        const disp = document.getElementById(b.disp);
        const val = parseFloat(el.value);
        SIM_PARAMS[b.param] = val;
        disp.textContent = b.fmt ? b.fmt(val) : `${val}${b.unit}`;
      });

      SIM_PARAMS.waterDensity = defaultPreset.waterDensity;

      if (this.isPaused) {
        this.isPaused = false;
        btnPause.textContent = "Pause";
        btnPause.classList.remove("btn-primary");
        btnPause.classList.add("btn-secondary");
      }

      this.onParamsChanged();
    });

    // 3. Settle to Equilibrium button
    document.getElementById("btn-reset-phys").addEventListener("click", () => {
      this.physics.resetToEquilibrium();
    });

    // 4. View Mode Toggle
    const btnViewToggle = document.getElementById("btn-view-toggle");
    btnViewToggle.addEventListener("click", () => {
      if (typeof THREE === "undefined") {
        alert("Three.js library is not loaded. The 2D sagittal cross-section mode provides full interactive physics and hydrostatic modeling.");
        return;
      }
      this.is3DMode = !this.is3DMode;
      if (this.is3DMode) {
        this.canvas2D.style.display = "none";
        this.container3D.style.display = "block";
        this.renderer3D.init();
        btnViewToggle.textContent = "Switch to 2D Cross-Section";
      } else {
        this.canvas2D.style.display = "block";
        this.container3D.style.display = "none";
        btnViewToggle.textContent = "Switch to 3D View";
      }
    });
  }

  onParamsChanged() {
    this.model.update(SIM_PARAMS);
    this.hydro = new HydrostaticSolver(this.model);
    this.physics.model = this.model;
    this.physics.hydro = this.hydro;
    this.physics.resetToEquilibrium();
    if (this.is3DMode && this.renderer3D.initialized) {
      this.renderer3D.buildShellMesh();
    }
    this.updateHUD();
  }

  bindMouseInteractions() {
    const canvas = this.canvas2D;
    const forceOverlay = document.getElementById("force-overlay");
    const forceVal = document.getElementById("force-val");
    const torqueVal = document.getElementById("torque-val");

    const getMousePos = (e) => {
      const rect = canvas.getBoundingClientRect();
      return {
        x: e.clientX - rect.left,
        y: e.clientY - rect.top
      };
    };

    canvas.addEventListener("mousedown", (e) => {
      const mouse = getMousePos(e);
      const wMouse = this.renderer2D.screenToWorld(mouse.x, mouse.y);
      const hitResult = this.renderer2D.hitTest(wMouse.x, wMouse.y);

      if (hitResult.hit) {
        this.physics.setDragStart(hitResult.bodyPt, wMouse);
        forceOverlay.classList.remove("hidden");
      }
    });

    window.addEventListener("mousemove", (e) => {
      if (this.physics.isDragging) {
        const mouse = getMousePos(e);
        const wMouse = this.renderer2D.screenToWorld(mouse.x, mouse.y);
        this.physics.setDragMove(wMouse);

        // Update force overlay values
        const fN = this.physics.appliedForce.mag;
        const fGf = fN * 101.97;
        forceVal.textContent = `${fN.toFixed(2)} N (${fGf.toFixed(0)} gf)`;
        torqueVal.textContent = `${this.physics.appliedTorque.toFixed(3)} N·m`;
      }
    });

    window.addEventListener("mouseup", () => {
      if (this.physics.isDragging) {
        this.physics.setDragEnd();
        forceOverlay.classList.add("hidden");
      }
    });

    // Touch events for mobile / touch devices
    canvas.addEventListener("touchstart", (e) => {
      if (e.touches.length === 1) {
        const t = e.touches[0];
        const rect = canvas.getBoundingClientRect();
        const wMouse = this.renderer2D.screenToWorld(t.clientX - rect.left, t.clientY - rect.top);
        const hitResult = this.renderer2D.hitTest(wMouse.x, wMouse.y);
        if (hitResult.hit) {
          e.preventDefault();
          this.physics.setDragStart(hitResult.bodyPt, wMouse);
          forceOverlay.classList.remove("hidden");
        }
      }
    }, { passive: false });

    window.addEventListener("touchmove", (e) => {
      if (this.physics.isDragging && e.touches.length > 0) {
        e.preventDefault();
        const t = e.touches[0];
        const rect = canvas.getBoundingClientRect();
        const wMouse = this.renderer2D.screenToWorld(t.clientX - rect.left, t.clientY - rect.top);
        this.physics.setDragMove(wMouse);

        const fN = this.physics.appliedForce.mag;
        const fGf = fN * 101.97;
        forceVal.textContent = `${fN.toFixed(2)} N (${fGf.toFixed(0)} gf)`;
        torqueVal.textContent = `${this.physics.appliedTorque.toFixed(3)} N·m`;
      }
    }, { passive: false });

    window.addEventListener("touchend", () => {
      if (this.physics.isDragging) {
        this.physics.setDragEnd();
        forceOverlay.classList.add("hidden");
      }
    });

    // Canvas zoom & pan
    canvas.addEventListener("wheel", (e) => {
      e.preventDefault();
      const zoomFactor = e.deltaY < 0 ? 1.08 : 0.92;
      this.renderer2D.zoom = Math.max(0.4, Math.min(3.0, this.renderer2D.zoom * zoomFactor));
    }, { passive: false });
  }

  updateHUD() {
    const eq = this.hydro.solveEquilibrium();
    const m = this.model;

    const hudStatus = document.getElementById("hud-status");
    const hudOrientation = document.getElementById("hud-orientation");
    const hudCom = document.getElementById("hud-com");
    const hudCob = document.getElementById("hud-cob");
    const hudStability = document.getElementById("hud-stability");
    const hudMoi = document.getElementById("hud-moi");
    const hudMass = document.getElementById("hud-mass");

    if (eq.isFloating) {
      hudStatus.textContent = "Positive (Floats)";
      hudStatus.className = "hud-value status-badge";
      hudOrientation.textContent = `${eq.angleDeg.toFixed(1)}°`;
      hudCob.textContent = `X: ${eq.cob.x.toFixed(0)}mm, Y: ${eq.cob.y.toFixed(0)}mm`;
      hudStability.textContent = `${eq.sepDistance.toFixed(1)}mm (St=${eq.stabilityIndex.toFixed(2)})`;
    } else {
      hudStatus.textContent = "Negative (Sinks)";
      hudStatus.className = "hud-value status-badge sinking";
      hudOrientation.textContent = "N/A (Benthic)";
      hudCob.textContent = "Submerged Base";
      hudStability.textContent = "Unstable (Sunk)";
    }

    hudCom.textContent = `X: ${m.com.x.toFixed(0)}mm, Y: ${m.com.y.toFixed(0)}mm`;
    hudMoi.textContent = `${m.Izz_SI.toFixed(4)} kg·m² (k=${m.k_gyr.toFixed(0)}mm)`;
    hudMass.textContent = `M: ${(m.totalMass).toFixed(0)}g | Shell: ${m.shellMass.toFixed(0)}g | Siph: ${m.siphMass.toFixed(0)}g`;
  }

  animate(currentTime) {
    const dt = Math.min(0.033, (currentTime - this.lastTime) * 0.001);
    this.lastTime = currentTime;

    // Advance physics only if not paused
    if (!this.isPaused) {
      this.physics.step(dt);
    }

    // Render active viewport
    if (this.is3DMode) {
      this.renderer3D.render();
    } else {
      this.renderer2D.render();
    }

    requestAnimationFrame(this.animate);
  }
}

// Initialize simulation on DOM loaded
window.addEventListener("DOMContentLoaded", () => {
  window.simApp = new SimulationApp();
});
