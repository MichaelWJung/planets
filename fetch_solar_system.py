#!/usr/bin/env python3
"""
fetch_solar_system.py

Queries the JPL Horizons REST API for every body in a curated list,
fetching state vectors (x, y, z, vx, vy, vz) and GM (gravitational
parameter = G * mass) at a single epoch, then writes a CSV file
suitable for loading directly into a C++ simulation.

Output columns:
  id, name, GM_km3_s2, x_km, y_km, z_km, vx_km_s, vy_km_s, vz_km_s

Coordinate frame : ECLIPJ2000, solar system barycenter origin (500@0)
Units            : km for position, km/s for velocity, km^3/s^2 for GM

Usage:
  pip install requests
  python fetch_solar_system.py                        # J2000 epoch
  python fetch_solar_system.py --epoch "2024-Jan-1"
  python fetch_solar_system.py --out myfile.csv

Notes:
  - GM=0.0 means Horizons did not report a GM value for that body.
    Masses for those can be looked up at:
      https://ssd.jpl.nasa.gov/sats/phys_par/
      https://ssd.jpl.nasa.gov/sb/phys_par.html
  - Many recently discovered small irregular moons of Jupiter/Saturn
    have provisional IDs and Horizons may lack ephemerides at J2000.
    The script reports failures; some are expected for tiny bodies.
  - The asteroid list is ordered roughly by diameter (largest first).
"""

import re
import csv
import time
import argparse
import requests
from datetime import datetime

# ===========================================================================
# BODY LIST
# Each entry: (horizons_id, human_readable_name)
#
# Major body integer IDs:
#   10        = Sun
#   X99       = planet center  (199=Mercury, 299=Venus, etc.)
#   X01..X99  = moons (501=Io, 601=Mimas, etc.)
#   9XX       = Pluto system
#
# Small bodies: "N;" format (e.g. "1;" = Ceres)
# ===========================================================================

BODIES = [

    # -------------------------------------------------------------------
    # SUN
    # -------------------------------------------------------------------
    (10,    "Sun"),

    # -------------------------------------------------------------------
    # PLANETS
    # -------------------------------------------------------------------
    (199,   "Mercury"),
    (299,   "Venus"),
    (399,   "Earth"),
    (499,   "Mars"),
    (599,   "Jupiter"),
    (699,   "Saturn"),
    (799,   "Uranus"),
    (899,   "Neptune"),

    # -------------------------------------------------------------------
    # EARTH MOON
    # -------------------------------------------------------------------
    (301,   "Moon"),

    # -------------------------------------------------------------------
    # MARS MOONS  (2)
    # -------------------------------------------------------------------
    (401,   "Phobos"),
    (402,   "Deimos"),

    # -------------------------------------------------------------------
    # JUPITER MOONS
    # 95 confirmed moons as of 2024. Horizons covers named + numbered
    # moons with IDs 501-572 and beyond.
    # -------------------------------------------------------------------
    (501,   "Io"),
    (502,   "Europa"),
    (503,   "Ganymede"),
    (504,   "Callisto"),
    (505,   "Amalthea"),
    (506,   "Himalia"),
    (507,   "Elara"),
    (508,   "Pasiphae"),
    (509,   "Sinope"),
    (510,   "Lysithea"),
    (511,   "Carme"),
    (512,   "Ananke"),
    (513,   "Leda"),
    (514,   "Thebe"),
    (515,   "Adrastea"),
    (516,   "Metis"),
    (517,   "Callirrhoe"),
    (518,   "Themisto"),
    (519,   "Megaclite"),
    (520,   "Taygete"),
    (521,   "Chaldene"),
    (522,   "Harpalyke"),
    (523,   "Kalyke"),
    (524,   "Iocaste"),
    (525,   "Erinome"),
    (526,   "Isonoe"),
    (527,   "Praxidike"),
    (528,   "Autonoe"),
    (529,   "Thyone"),
    (530,   "Hermippe"),
    (531,   "Aitne"),
    (532,   "Eurydome"),
    (533,   "Euanthe"),
    (534,   "Euporie"),
    (535,   "Orthosie"),
    (536,   "Sponde"),
    (537,   "Kale"),
    (538,   "Pasithee"),
    (539,   "Hegemone"),
    (540,   "Mneme"),
    (541,   "Aoede"),
    (542,   "Thelxinoe"),
    (543,   "Arche"),
    (544,   "Kallichore"),
    (545,   "Helike"),
    (546,   "Carpo"),
    (547,   "Eukelade"),
    (548,   "Cyllene"),
    (549,   "Kore"),
    (550,   "Herse"),
    (551,   "Dia"),
    (552,   "Pandia"),
    (553,   "Ersa"),
    (554,   "Philophrosyne"),
    (555,   "Eupheme"),
    (556,   "Valetudo"),
    (557,   "Eirene"),
    (558,   "Philophrosyne-2"),
    (560,   "S/2010 J1"),
    (561,   "S/2010 J2"),
    (562,   "S/2011 J1"),
    (563,   "S/2011 J2"),
    (564,   "S/2016 J1"),
    (565,   "S/2017 J1"),
    (566,   "S/2017 J2"),
    (567,   "S/2017 J3"),
    (568,   "S/2017 J5"),
    (569,   "S/2017 J6"),
    (570,   "S/2017 J7"),
    (571,   "S/2017 J8"),
    (572,   "S/2017 J9"),

    # -------------------------------------------------------------------
    # SATURN MOONS
    # 146 confirmed moons as of 2024. Named IDs 601-668+.
    # -------------------------------------------------------------------
    (601,   "Mimas"),
    (602,   "Enceladus"),
    (603,   "Tethys"),
    (604,   "Dione"),
    (605,   "Rhea"),
    (606,   "Titan"),
    (607,   "Hyperion"),
    (608,   "Iapetus"),
    (609,   "Phoebe"),
    (610,   "Janus"),
    (611,   "Epimetheus"),
    (612,   "Helene"),
    (613,   "Telesto"),
    (614,   "Calypso"),
    (615,   "Atlas"),
    (616,   "Prometheus"),
    (617,   "Pandora"),
    (618,   "Pan"),
    (619,   "Ymir"),
    (620,   "Paaliaq"),
    (621,   "Tarvos"),
    (622,   "Ijiraq"),
    (623,   "Suttungr"),
    (624,   "Kiviuq"),
    (625,   "Mundilfari"),
    (626,   "Albiorix"),
    (627,   "Skathi"),
    (628,   "Erriapus"),
    (629,   "Siarnaq"),
    (630,   "Thrymr"),
    (631,   "Narvi"),
    (632,   "Methone"),
    (633,   "Pallene"),
    (634,   "Polydeuces"),
    (635,   "Daphnis"),
    (636,   "Aegaeon"),
    (637,   "Aegir"),
    (638,   "Bebhionn"),
    (639,   "Bergelmir"),
    (640,   "Bestla"),
    (641,   "Farbauti"),
    (642,   "Fenrir"),
    (643,   "Fornjot"),
    (644,   "Hati"),
    (645,   "Hyrrokkin"),
    (646,   "Kari"),
    (647,   "Loge"),
    (648,   "Skoll"),
    (649,   "Anthe"),
    (650,   "Surtur"),
    (651,   "Greip"),
    (652,   "Jarnsaxa"),
    (653,   "Tarqeq"),
    (655,   "S/2004 S7"),
    (656,   "S/2004 S12"),
    (657,   "S/2004 S13"),
    (658,   "S/2004 S17"),
    (659,   "S/2006 S1"),
    (660,   "S/2006 S3"),
    (661,   "S/2007 S2"),
    (662,   "S/2007 S3"),
    (663,   "S/2009 S1"),
    (664,   "Gridr"),
    (665,   "Angrboda"),
    (666,   "Skrymir"),
    (667,   "Gerd"),
    (668,   "S/2004 S26"),

    # -------------------------------------------------------------------
    # URANUS MOONS  (28 named; IDs 701-728)
    # -------------------------------------------------------------------
    (701,   "Ariel"),
    (702,   "Umbriel"),
    (703,   "Titania"),
    (704,   "Oberon"),
    (705,   "Miranda"),
    (706,   "Cordelia"),
    (707,   "Ophelia"),
    (708,   "Bianca"),
    (709,   "Cressida"),
    (710,   "Desdemona"),
    (711,   "Juliet"),
    (712,   "Portia"),
    (713,   "Rosalind"),
    (714,   "Belinda"),
    (715,   "Puck"),
    (716,   "Caliban"),
    (717,   "Sycorax"),
    (718,   "Prospero"),
    (719,   "Setebos"),
    (720,   "Stephano"),
    (721,   "Trinculo"),
    (722,   "Francisco"),
    (723,   "Margaret"),
    (724,   "Ferdinand"),
    (725,   "Perdita"),
    (726,   "Mab"),
    (727,   "Cupid"),
    (728,   "S/2023 U1"),

    # -------------------------------------------------------------------
    # NEPTUNE MOONS  (16 named; IDs 801-816)
    # -------------------------------------------------------------------
    (801,   "Triton"),
    (802,   "Nereid"),
    (803,   "Naiad"),
    (804,   "Thalassa"),
    (805,   "Despina"),
    (806,   "Galatea"),
    (807,   "Larissa"),
    (808,   "Proteus"),
    (809,   "Halimede"),
    (810,   "Sao"),
    (811,   "Laomedeia"),
    (812,   "Psamathe"),
    (813,   "Neso"),
    (814,   "Hippocamp"),
    (815,   "S/2021 N1"),
    (816,   "S/2002 N4"),

    # -------------------------------------------------------------------
    # PLUTO SYSTEM  (5 moons)
    # -------------------------------------------------------------------
    (999,   "Pluto"),
    (901,   "Charon"),
    (902,   "Nix"),
    (903,   "Hydra"),
    (904,   "Kerberos"),
    (905,   "Styx"),

    # -------------------------------------------------------------------
    # DWARF PLANETS
    # IAU-official (5): Ceres (below in asteroids), Pluto (above),
    # Eris, Haumea, Makemake, Quaoar (2024)
    # Strong candidates: Sedna, Orcus, Gonggong, Varuna, Ixion, etc.
    # -------------------------------------------------------------------
    ("136199;", "Eris"),
    ("136108;", "Haumea"),
    ("136472;", "Makemake"),
    ("50000;",  "Quaoar"),
    ("90377;",  "Sedna"),
    ("90482;",  "Orcus"),
    ("225088;", "Gonggong"),
    ("20000;",  "Varuna"),
    ("28978;",  "Ixion"),
    ("55565;",  "2002 AW197"),
    ("55637;",  "2002 UX25"),
    ("84522;",  "2002 TC302"),
    ("120347;", "Salacia"),
    ("174567;", "Varda"),
    ("523794;", "2015 RR245"),
    ("471143;", "2010 EK139"),  # Dziewanna
    ("444030;", "2004 UX10"),
    ("308193;", "2005 QU182"),

    # Moons of dwarf planets (Horizons major body IDs where available)
    (120348,  "Hi'iaka"),       # Haumea I
    (120349,  "Namaka"),        # Haumea II
    (90483,   "Vanth"),         # Orcus I
    (50001,   "Weywot"),        # Quaoar I

    # -------------------------------------------------------------------
    # ASTEROIDS — Top ~100+ by size/mass
    # Ordered roughly largest-first by diameter.
    # All use "N;" small-body format.
    # -------------------------------------------------------------------
    # The big four (contain >50% of belt mass)
    ("1;",    "Ceres"),
    ("4;",    "Vesta"),
    ("2;",    "Pallas"),
    ("10;",   "Hygiea"),
    # 5-30 by diameter (>200 km)
    ("704;",  "Interamnia"),
    ("52;",   "Europa-ast"),
    ("511;",  "Davida"),
    ("87;",   "Sylvia"),
    ("65;",   "Cybele"),
    ("15;",   "Eunomia"),
    ("16;",   "Psyche"),
    ("31;",   "Euphrosyne"),
    ("107;",  "Camilla"),
    ("324;",  "Bamberga"),
    ("3;",    "Juno"),
    ("88;",   "Thisbe"),
    ("451;",  "Patientia"),
    ("532;",  "Herculina"),
    ("48;",   "Doris"),
    ("375;",  "Ursula"),
    ("409;",  "Aspasia"),
    ("423;",  "Diotima"),
    ("24;",   "Themis"),
    ("128;",  "Nemesis"),
    ("45;",   "Eugenia"),
    ("121;",  "Hermione"),
    ("130;",  "Elektra"),
    # 31-60 by size (~150-200 km)
    ("19;",   "Fortuna"),
    ("41;",   "Daphne"),
    ("7;",    "Iris"),
    ("6;",    "Hebe"),
    ("9;",    "Metis"),
    ("22;",   "Kalliope"),
    ("29;",   "Amphitrite"),
    ("43;",   "Ariadne"),
    ("46;",   "Hestia"),
    ("47;",   "Aglaja"),
    ("50;",   "Virginia"),
    ("54;",   "Alexandra"),
    ("59;",   "Elpis"),
    ("62;",   "Erato"),
    ("70;",   "Panopaea"),
    ("76;",   "Freia"),
    ("78;",   "Diana"),
    ("85;",   "Io-ast"),
    ("93;",   "Minerva"),
    ("94;",   "Aurora"),
    ("95;",   "Arethusa"),
    ("96;",   "Aegle"),
    ("104;",  "Klymene"),
    ("106;",  "Dione-ast"),
    ("109;",  "Felicitas"),
    ("117;",  "Lomia"),
    ("120;",  "Lachesis"),
    ("129;",  "Antigone"),
    ("139;",  "Juewa"),
    ("140;",  "Siwa"),
    # 61-100 by size (~120-150 km)
    ("145;",  "Adeona"),
    ("150;",  "Nuwa"),
    ("154;",  "Bertha"),
    ("159;",  "Aemilia"),
    ("162;",  "Laurentia"),
    ("168;",  "Sibylla"),
    ("173;",  "Ino"),
    ("175;",  "Andromache"),
    ("181;",  "Eucharis"),
    ("185;",  "Eunike"),
    ("187;",  "Lamberta"),
    ("194;",  "Prokne"),
    ("196;",  "Philomela"),
    ("200;",  "Dynamene"),
    ("211;",  "Isolda"),
    ("212;",  "Medea"),
    ("216;",  "Kleopatra"),
    ("221;",  "Eos"),
    ("225;",  "Henrietta"),
    ("233;",  "Asterope"),
    ("238;",  "Hypatia"),
    ("241;",  "Germania"),
    ("247;",  "Eukrate"),
    ("250;",  "Bettina"),
    ("259;",  "Aletheia"),
    ("283;",  "Emma"),
    ("308;",  "Polyxo"),
    ("334;",  "Chicago"),
    ("344;",  "Desiderata"),
    ("356;",  "Liguria"),
    ("372;",  "Palma"),
    ("381;",  "Myrrha"),
    ("386;",  "Siegena"),
    ("388;",  "Charybdis"),
    ("407;",  "Arachne"),
    ("410;",  "Chloris"),
    ("419;",  "Aurelia"),
    ("420;",  "Bertholda"),
    ("444;",  "Gyptis"),
    ("449;",  "Hamburga"),
    ("476;",  "Hedwig"),
    ("488;",  "Kreusa"),
    ("489;",  "Comacina"),
    # Notable / famous objects
    ("433;",  "Eros"),
    ("243;",  "Ida"),
    ("253;",  "Mathilde"),
    ("99942;","Apophis"),
    ("101955;","Bennu"),
    ("162173;","Ryugu"),
    ("25143;","Itokawa"),
    ("216;",  "Kleopatra"),   # already listed; will be deduped
    # Jupiter Trojans (large ones)
    ("624;",  "Hektor"),
    ("617;",  "Patroclus"),
    ("588;",  "Achilles"),
    ("659;",  "Nestor"),
    ("1172;", "Aneas"),
    ("1143;", "Odysseus"),
    ("911;",  "Agamemnon"),
    ("1437;", "Diomedes"),
    ("884;",  "Priamus"),
    # Centaurs
    ("2060;", "Chiron"),
    ("5145;", "Pholus"),
    ("10199;","Chariklo"),
    ("7066;", "Nessus"),
    ("54598;","Bienor"),
    ("52872;","Okyrhoe"),
]

# De-duplicate while preserving order
seen = set()
BODIES_DEDUPED = []
for item in BODIES:
    key = str(item[0])
    if key not in seen:
        seen.add(key)
        BODIES_DEDUPED.append(item)
BODIES = BODIES_DEDUPED


# ===========================================================================
# Horizons API helpers
# ===========================================================================
API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

def _epoch_to_jd(epoch: str) -> float:
    """Convert a Horizons-style epoch string to a Julian Day Number.

    Accepts: 'YYYY-Mmm-D HH:MM' or 'YYYY-Mmm-D'.
    """
    for fmt in ("%Y-%b-%d %H:%M", "%Y-%b-%d"):
        try:
            dt = datetime.strptime(epoch, fmt)
            break
        except ValueError:
            pass
    else:
        raise ValueError(f"Unrecognised epoch format: {epoch!r}")

    y, m, d = dt.year, dt.month, dt.day
    frac = (dt.hour * 60 + dt.minute) / 1440.0
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 4
    jd = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5 + frac
    return jd

def build_params(command: str, epoch: str) -> dict:
    jd_start = _epoch_to_jd(epoch)
    jd_stop  = jd_start + 1.0
    return {
        "format":      "text",
        "COMMAND":     command,
        "OBJ_DATA":    "YES",
        "MAKE_EPHEM":  "YES",
        "EPHEM_TYPE":  "VECTORS",
        "CENTER":      "500@0",
        "REF_PLANE":   "ECLIPTIC",
        "REF_SYSTEM":  "J2000",
        "START_TIME":  f"JD{jd_start:.6f}",
        "STOP_TIME":   f"JD{jd_stop:.6f}",
        "STEP_SIZE":   "1d",
        "VEC_TABLE":   "2",
        "VEC_CORR":    "NONE",
        "OUT_UNITS":   "KM-S",
        "CSV_FORMAT":  "NO",
        "VEC_DELTA_T": "NO",
    }

def parse_gm(text: str) -> float:
    patterns = [
        r"GM\s*\(?\s*km\^3\s*/\s*s\^2\s*\)?\s*[=,]\s*([\d.eE+\-]+)",
        r"GM\s*=\s*([\d.eE+\-]+)",
        r"GM,\s*km\^3/s\^2\s*[=:]\s*([\d.eE+\-]+)",
    ]
    for pat in patterns:
        m = re.search(pat, text, re.IGNORECASE)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                pass
    return 0.0

def parse_vectors(text: str):
    soe = text.find("$$SOE")
    eoe = text.find("$$EOE")
    if soe == -1 or eoe == -1:
        return None
    block = text[soe:eoe]

    def grab(label):
        m = re.search(rf"{label}\s*=\s*([\-\d.eE+]+)", block, re.IGNORECASE)
        return float(m.group(1)) if m else None

    x  = grab("X ")
    y  = grab("Y ")
    z  = grab("Z ")
    vx = grab("VX")
    vy = grab("VY")
    vz = grab("VZ")

    if None in (x, y, z, vx, vy, vz):
        return None
    return x, y, z, vx, vy, vz

def query_body(command, epoch: str, retries: int = 3):
    params = build_params(str(command), epoch)
    for attempt in range(retries):
        try:
            r = requests.get(API_URL, params=params, timeout=30)
            r.raise_for_status()
            text = r.text

            bad = ("No ephemeris for target", "cannot be",
                   "not found", "multiple major-bodies",
                   "ambiguous", "No matches found")
            if any(b.lower() in text.lower() for b in bad):
                return None

            gm  = parse_gm(text)
            vec = parse_vectors(text)
            if vec is None:
                return None
            return (gm,) + vec

        except requests.exceptions.RequestException as e:
            print(f"\n    [error] attempt {attempt+1}/{retries}: {e}")
            time.sleep(2 ** attempt)
    return None


# ===========================================================================
# Main
# ===========================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Fetch solar system state vectors from JPL Horizons")
    parser.add_argument("--epoch", default="2000-Jan-1 12:00",
        help="Epoch in Horizons format (default: 2000-Jan-1 12:00 = J2000)")
    parser.add_argument("--out", default="solar_system.csv",
        help="Output CSV filename (default: solar_system.csv)")
    parser.add_argument("--delay", type=float, default=0.6,
        help="Seconds between API calls (default: 0.6)")
    args = parser.parse_args()

    total = len(BODIES)
    print(f"Epoch   : {args.epoch}")
    print(f"Output  : {args.out}")
    print(f"Targets : {total} unique bodies queued")
    print(f"Est. time: ~{total * args.delay / 60:.0f} minutes")
    print()

    rows   = []
    failed = []

    for idx, (body_id, name) in enumerate(BODIES, 1):
        label = f"[{idx:3d}/{total}] {name:32s} (id={body_id})"
        print(label, end=" ... ", flush=True)
        result = query_body(body_id, args.epoch)
        if result is None:
            print("SKIP")
            failed.append((body_id, name))
        else:
            gm, x, y, z, vx, vy, vz = result
            gm_str = f"GM={gm:.4g}" if gm else "GM=?"
            print(f"OK  {gm_str}")
            rows.append({
                "id":        str(body_id),
                "name":      name,
                "GM_km3_s2": gm,
                "x_km":      x,
                "y_km":      y,
                "z_km":      z,
                "vx_km_s":   vx,
                "vy_km_s":   vy,
                "vz_km_s":   vz,
            })
        time.sleep(args.delay)

    # Write output
    fields = ["id", "name", "GM_km3_s2",
              "x_km", "y_km", "z_km",
              "vx_km_s", "vy_km_s", "vz_km_s"]

    with open(args.out, "w", newline="") as f:
        f.write(f"# JPL Horizons solar system state vectors\n")
        f.write(f"# Epoch    : {args.epoch}\n")
        f.write(f"# Frame    : ECLIPJ2000, solar system barycenter (500@0)\n")
        f.write(f"# Position : km\n")
        f.write(f"# Velocity : km/s\n")
        f.write(f"# GM       : km^3/s^2  (GM=0 = not reported by Horizons)\n")
        f.write(f"# Fetched  : {len(rows)} OK, {len(failed)} skipped\n")
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print()
    print("=" * 65)
    print(f"Wrote {len(rows)} bodies to '{args.out}'")
    print(f"Skipped {len(failed)} (no Horizons ephemeris at this epoch)")

    if failed:
        print(f"\nSkipped bodies:")
        for body_id, name in failed:
            print(f"  {name:32s}  id={body_id}")

    no_gm = [r for r in rows if r["GM_km3_s2"] == 0.0]
    if no_gm:
        print(f"\n{len(no_gm)} bodies have GM=0. Look up masses at:")
        print(f"  https://ssd.jpl.nasa.gov/sats/phys_par/    (moons)")
        print(f"  https://ssd.jpl.nasa.gov/sb/phys_par.html  (small bodies)")

if __name__ == "__main__":
    main()
