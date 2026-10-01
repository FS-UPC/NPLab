# Experimental dataset GOx / M10pH7 -- three measurement sessions, three separate ZIPs

This material corresponds to Fig. 2 of the article. This folder contains **three independent ZIP files**, one per measurement session: `F61.zip`,
`F62.zip`, and `F63D.zip`.

NPLab loads the calibration together with the rest of a session's files, and that
calibration stays active (set) to process anything loaded afterwards. Since F61, F62,
and F63D are three independent measurement sessions, each has its **own calibration
curve**, fitted from its own CCN standards -- these curves are not equivalent to one
another (different slope and intercept). If files from more than one session are loaded
together in NPLab, the software will process one session's spectra using another
session's calibration, producing incorrect concentration and adsorption-capacity values.

For this reason, each session must be loaded **separately** into NPLab -- one ZIP at a
time, each with its own calibration active -- never mixing files from two different ZIPs
within the same processing session.

## What each session represents

- **F61** -- kinetics at medium/long contact times (t = 5-300 min), characterising the
  adsorption equilibrium plateau.
- **F62** -- kinetics at the shortest contact times (t = 1, 2 min), designed specifically
  to determine how quickly equilibrium is reached. Also includes the Control (CTRLN)
  files used to apply Co-correction (against the C000 control).
- **F63D** -- the full isotherm (N = 5-400 mg/L, 9 levels), with its own calibration,
  blanks, drift, and Co-correction control files.

## Contents of each ZIP

- `F61.zip` (31 files): calibration (CCN x8), blanks (BLK x2), drift (DRF x7), kinetics
  (AKN x14, t=5,15,60,120,180,240,300 min), and `kinetics_data_summary_F61.csv` with the
  expected (t, q, q_std) values.
- `F62.zip` (75 files): calibration (CCN x15), blanks (BLK x30), Co-correction controls
  (CTRLN x9), drift (DRF x12), kinetics (AKN x6, t=1,2 min), and
  `kinetics_data_summary_F62.csv`.
- `F63D.zip` (109 files): calibration (CCN x23), blanks (BLK x26), Co-correction controls
  (CTRLN x21), drift (DRF x12), isotherm (AIN x27, N=5-400 mg/L), and
  `isotherm_data_summary_F63D.csv`.

## How to use in NPLab

1. Open NPLab and load **one of the three ZIPs** (e.g., `F61.zip`).
2. Go to the **Calibration** tab and load/confirm that ZIP's `CCN*.txt` files -- this
   sets that session's calibration.
3. If the ZIP includes `CTRLN*.txt` files (F62 and F63D), go to the **Control** tab,
   load them, and enable Co-correction against the C000 control.
4. Go to the corresponding tab (**Kinetics** for F61/F62, **Isotherms** for F63D) and
   load that same ZIP's `AKN*.txt` or `AIN*.txt` files.
5. Check the results obtained against the reference CSV included in that ZIP.
6. Close/clear the session before repeating the process with the next ZIP.

Reprocessing each session separately should reproduce, within NPLab's own reported
precision, the numbers quoted in the article: kinetics plateau 51.1+-0.9 mg/g (CV=1.8%),
PFO/PSO R^2~=0 (Fig. 2a); Langmuir qmax=349.1 mg/g, KL=0.00216 L/mg, R^2=0.989,
Freundlich 1/n=0.799, R^2=0.983 (Fig. 2b).
