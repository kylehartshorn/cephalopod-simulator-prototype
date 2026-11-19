/**
 * Paleontological Morphometric Presets for Ectocochleate Cephalopods
 * References:
 * - Trueman (1941) Q. J. Geol. Soc. Lond.
 * - Westermann (1977) Paleobiology
 * - Saunders & Shapiro (1986) Paleobiology
 * - Peterman et al. (2020) Scientific Reports
 */

const CEPHALOPOD_PRESETS = {
  orthoceras: {
    name: "Orthoceras regulare",
    age: "Middle Ordovician – Silurian",
    shellLength: 600,
    curvature: 0.0,
    apertureDiam: 80,
    livingChamberLength: 180,
    shellThickness: 2.0,
    septaThickness: 1.0,
    chamberDist: 10,
    siphuncleDiam: 10,
    siphunclePos: 50, // Central
    siphuncleDensity: 2.70,
    shellDensity: 2.70,
    chamberGasDensity: 0.0012, // Standard air ~1.2 kg/m3 (0.0012 g/cm3)
    weightFactor: 1.00,
    cameralFluid: 0,
    waterDensity: 1.025
  },
  endoceras: {
    name: "Endoceras proteiforme",
    age: "Middle Ordovician",
    shellLength: 1500,
    curvature: 0.0,
    apertureDiam: 140,
    livingChamberLength: 450,
    shellThickness: 3.5,
    septaThickness: 1.8,
    chamberDist: 18,
    siphuncleDiam: 45, // Massive ventral siphuncle with calcified endocones
    siphunclePos: 0,   // Ventral margin
    siphuncleDensity: 2.70,
    shellDensity: 2.70,
    chamberGasDensity: 0.0012,
    weightFactor: 1.10, // Massive visceral mass
    cameralFluid: 0,
    waterDensity: 1.025
  },
  cyrtoceras: {
    name: "Cyrtoceras depressum",
    age: "Silurian – Devonian",
    shellLength: 420,
    curvature: 0.22,
    apertureDiam: 65,
    livingChamberLength: 135,
    shellThickness: 1.8,
    septaThickness: 0.9,
    chamberDist: 9,
    siphuncleDiam: 8,
    siphunclePos: 15, // Subventral
    siphuncleDensity: 2.70,
    shellDensity: 2.70,
    chamberGasDensity: 0.0012,
    weightFactor: 0.95,
    cameralFluid: 0,
    waterDensity: 1.025
  },
  gyronautilus: {
    name: "Gyronautilus",
    age: "Devonian – Carboniferous",
    shellLength: 550,
    curvature: 0.50,
    apertureDiam: 75,
    livingChamberLength: 170,
    shellThickness: 2.0,
    septaThickness: 1.0,
    chamberDist: 11,
    siphuncleDiam: 10,
    siphunclePos: 40, // Subcentral
    siphuncleDensity: 2.70,
    shellDensity: 2.70,
    chamberGasDensity: 0.0012,
    weightFactor: 1.00,
    cameralFluid: 0,
    waterDensity: 1.025
  },
  nautilus: {
    name: "Nautilus pompilius",
    age: "Extant (Relict Nautilida)",
    shellLength: 680,
    curvature: 0.95,
    apertureDiam: 90,
    livingChamberLength: 220,
    shellThickness: 1.6,
    septaThickness: 0.8,
    chamberDist: 12,
    siphuncleDiam: 7,
    siphunclePos: 45, // Subcentral
    siphuncleDensity: 2.93,
    shellDensity: 2.93, // Aragonite conch
    chamberGasDensity: 0.0012,
    weightFactor: 1.02,
    cameralFluid: 5,    // Pellicle and siphuncular fluid
    waterDensity: 1.025
  },
  baculites: {
    name: "Baculites compressus",
    age: "Late Cretaceous",
    shellLength: 800,
    curvature: 0.03,
    apertureDiam: 55,
    livingChamberLength: 260,
    shellThickness: 1.3,
    septaThickness: 0.5,
    chamberDist: 12,
    siphuncleDiam: 5,
    siphunclePos: 0, // Marginal ventral
    siphuncleDensity: 2.93,
    shellDensity: 2.93, // Aragonite
    chamberGasDensity: 0.0012,
    weightFactor: 1.00,
    cameralFluid: 0,
    waterDensity: 1.025
  },
  dactylioceras: {
    name: "Dactylioceras commune",
    age: "Early Jurassic",
    shellLength: 520,
    curvature: 0.85,
    apertureDiam: 48,
    livingChamberLength: 195,
    shellThickness: 1.2,
    septaThickness: 0.5,
    chamberDist: 8,
    siphuncleDiam: 5,
    siphunclePos: 0, // Marginal ventral
    siphuncleDensity: 2.93,
    shellDensity: 2.93,
    chamberGasDensity: 0.0012,
    weightFactor: 1.04,
    cameralFluid: 0,
    waterDensity: 1.025
  }
};
