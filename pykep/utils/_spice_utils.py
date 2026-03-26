import spiceypy as pyspice
import pykep as pk
from pathlib import Path

def spice_version():
    """Retrieves the installed version of the NAIF spice toolkit.

    Returns:
        :class:`string`: SPICE toolbox version installed in the system.
    """
    return pyspice.tkvrsn("TOOLKIT")

def load_spice_kernels(paths):
    """Loads one or more kernels for use with the NAIF spice toolkit.

    Args:
        paths (:class:`string` or :class:`list`): The path (or paths) of the kernels to load
    """
    pyspice.furnsh(paths)

def unload_spice_kernels(paths):
    """Unloads from memory one or more kernels for use with the NAIF spice toolkit.

    Args:
        paths (:class:`string` or :class:`list`): The path (or paths) of the kernels to load
    """
    pyspice.unload(paths)

def inspect_spice_kernel(path):
    """Detects the objects contained in a specific kernel.
    This is useful to understand the possibile queries to make for a certain kernel.

    Args:
        path (:class:`string`): The path (or paths) of the kernels to load

    Returns:
        :class:`list`: the NAIF ids of the objects found in the kernels.
    """
    return list(pyspice.spkobj(path))

def extract_coverage_window(path, naifid, window_n = 0):
    """Extracts the *n*-th SPK coverage window for a given NAIF ID and returns it as PyKEP epochs.

    The SPK coverage start/end times returned by SPICE are ET (TDB) seconds past J2000. Here we map them to a **uniform**
    seconds-from-2000-01-01T00:00:00 (UTC-labeled origin) representation, so that the resulting PyKEP epoch is strictly
    1:1 with ET everywhere (this is not “UTC elapsed seconds” with leap seconds).

    The offset constant 43135.816087188054 [s] is ET0 = str2et("2000-01-01 00:00:00 UTC"), i.e. the ET value
    of the UTC-labeled origin of pykep mjd2000; thus mjd2000_days = (et - ET0) / 86400, implemented
    as (et + 43135.816087188054) * SEC2DAY.

    Args:
        path (:class:`string`): The path (or paths) of the kernel where the ephemerides are defined in windows.
        naifid (:class:`int`): The NAIF id.
        window_n (:class:`int`): The window to consider (in case ephs are defined in multiple windows of validity)
    
    Returns:
        :class:`pk.epoch`, :class:`pk.epoch`: Starting and final epoch of the *n*-th window
    """
    cover = pyspice.spkcov(path, naifid)
    start, end = pyspice.wnfetd(cover, window_n)
    start_mjd2000 = (start + 43135.816087188054) * pk.SEC2DAY # magic number is str2et("2000-01-01 00:00:00 UTC")
    end_mjd2000 = (end + 43135.816087188054) * pk.SEC2DAY # magic number is str2et("2000-01-01 00:00:00 UTC")
    #(we add/remove 1e-6 seconds as spice has better precision on the time representation and we want to avoid to go out of bunds)
    return pk.epoch(start_mjd2000)+1e-10, pk.epoch(end_mjd2000)-1e-10 

def epoch2utc(pykep_epoch):
    """Converts a PyKEP epoch (mjd2000) into a SPICE-formatted UTC calendar string.

    The input epoch is assumed to be a uniform time variable expressed as days since the UTC-labeled origin
    2000-01-01 00:00:00 (i.e. strictly 1:1 with SPICE ET everywhere; it is not “UTC elapsed seconds” with leap seconds).

    The constant 43135.816087188054 [s] is ET0 = str2et("2000-01-01 00:00:00 UTC"). We first reconstruct the ET value
    via et = mjd2000*86400 - ET0, then use SPICE et2utc to obtain the UTC calendar representation. A leapseconds kernel
    is loaded to enable correct ET↔UTC conversion.

    Args:
        pykep_epoch (:class:`pk.epoch`): PyKEP epoch whose .mjd2000 is days since 2000-01-01 00:00:00 (UTC-labeled origin).

    Returns:
        :class:`string`: UTC calendar string corresponding to the input epoch.
    """
    # We load the leapsecond kernel to allow exact conversions between 
    pk_path = Path(pk.__path__[0])
    kernel_leap = str(pk_path / "data" / "naif0012.tls")
    pk.utils.load_spice_kernels(kernel_leap)
    # We convert mjd2000 to et
    et = -43135.816087188054 + pykep_epoch.mjd2000 * pk.DAY2SEC
    return pyspice.et2utc(et, "C", 6)

def name2naifid(name):
    """Retreives the NAIF id of some planet/comet/spacecraft 

    Args:
        name (:class:`string`): The name of the planet/comet/spacecraft to be found

    Returns:
        :class:`int` the NAIF id of the the planet/comet/spacecraft
    """
    retval = pyspice.bodn2c(name)
    return retval

def framename2naifid(name):
    """Retreives the NAIF id of some reference frame

    Args:
        name (:class:`string`): The name of the reference frame to be found

    Returns:
        :class:`int` the NAIF id of the reference frame
    """
    retval = pyspice.irfnum(name)
    if retval == 0:
        raise NameError(name + " not found in the NAIF frames")
    return retval

def rotation_matrix(origin, destination, ep = pk.epoch(0)):
    """Rotation matrix between frames at epoch.

    Args:
        origin (:class:`string`): The name of the origin reference frame.

        destination (:class:`string`): The name of the destination reference frame.

        ep (:class:`~pykep.epoch`): Epoch. Defaults to pk.epoch(0).

    Returns:
        :class:`npumpy.ndarray`: The rotation matrix.
    """
    return pyspice.pxform(origin, destination, (ep-0.5).mjd2000*pk.DAY2SEC)


# These are taken from:
# https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html#NAIF%20Integer%20ID%20codes (13/10/2023)
naif_dict = dict(
    {
        0: "SOLAR SYSTEM BARYCENTER",
        1: "MERCURY BARYCENTER",
        2: "VENUS BARYCENTER",
        3: "EARTH BARYCENTER",
        4: "MARS BARYCENTER",
        5: "JUPITER BARYCENTER",
        6: "SATURN BARYCENTER",
        7: "URANUS BARYCENTER",
        8: "NEPTUNE BARYCENTER",
        9: "PLUTO BARYCENTER",
        10: "SUN",
        199: "MERCURY",
        299: "VENUS",
        399: "EARTH",
        301: "MOON",
        499: "MARS",
        401: "PHOBOS",
        402: "DEIMOS",
        599: "JUPITER",
        501: "IO",
        502: "EUROPA",
        503: "GANYMEDE",
        504: "CALLISTO",
        505: "AMALTHEA",
        506: "HIMALIA",
        507: "ELARA",
        508: "PASIPHAE",
        509: "SINOPE",
        510: "LYSITHEA",
        511: "CARME",
        512: "ANANKE",
        513: "LEDA",
        514: "THEBE",
        515: "ADRASTEA",
        516: "METIS",
        517: "CALLIRRHOE",
        518: "THEMISTO",
        519: "MEGACLITE",
        520: "TAYGETE",
        521: "CHALDENE",
        522: "HARPALYKE",
        523: "KALYKE",
        524: "IOCASTE",
        525: "ERINOME",
        526: "ISONOE",
        527: "PRAXIDIKE",
        528: "AUTONOE",
        529: "THYONE",
        530: "HERMIPPE",
        531: "AITNE",
        532: "EURYDOME",
        533: "EUANTHE",
        534: "EUPORIE",
        535: "ORTHOSIE",
        536: "SPONDE",
        537: "KALE",
        538: "PASITHEE",
        539: "HEGEMONE",
        540: "MNEME",
        541: "AOEDE",
        542: "THELXINOE",
        543: "ARCHE",
        544: "KALLICHORE",
        545: "HELIKE",
        546: "CARPO",
        547: "EUKELADE",
        548: "CYLLENE",
        549: "KORE",
        550: "HERSE",
        553: "DIA",
        699: "SATURN",
        601: "MIMAS",
        602: "ENCELADUS",
        603: "TETHYS",
        604: "DIONE",
        605: "RHEA",
        606: "TITAN",
        607: "HYPERION",
        608: "IAPETUS",
        609: "PHOEBE",
        610: "JANUS",
        611: "EPIMETHEUS",
        612: "HELENE",
        613: "TELESTO",
        614: "CALYPSO",
        615: "ATLAS",
        616: "PROMETHEUS",
        617: "PANDORA",
        618: "PAN",
        619: "YMIR",
        620: "PAALIAQ",
        621: "TARVOS",
        622: "IJIRAQ",
        623: "SUTTUNGR",
        624: "KIVIUQ",
        625: "MUNDILFARI",
        626: "ALBIORIX",
        627: "SKATHI",
        628: "ERRIAPUS",
        629: "SIARNAQ",
        630: "THRYMR",
        631: "NARVI",
        632: "METHONE",
        633: "PALLENE",
        634: "POLYDEUCES",
        635: "DAPHNIS",
        636: "AEGIR",
        637: "BEBHIONN",
        638: "BERGELMIR",
        639: "BESTLA",
        640: "FARBAUTI",
        641: "FENRIR",
        642: "FORNJOT",
        643: "HATI",
        644: "HYRROKKIN",
        645: "KARI",
        646: "LOGE",
        647: "SKOLL",
        648: "SURTUR",
        649: "ANTHE",
        650: "JARNSAXA",
        651: "GREIP",
        652: "TARQEQ",
        653: "AEGAEON",
        799: "URANUS",
        701: "ARIEL",
        702: "UMBRIEL",
        703: "TITANIA",
        704: "OBERON",
        705: "MIRANDA",
        706: "CORDELIA",
        707: "OPHELIA",
        708: "BIANCA",
        709: "CRESSIDA",
        710: "DESDEMONA",
        711: "JULIET",
        712: "PORTIA",
        713: "ROSALIND",
        714: "BELINDA",
        715: "PUCK",
        716: "CALIBAN",
        717: "SYCORAX",
        718: "PROSPERO",
        719: "SETEBOS",
        720: "STEPHANO",
        721: "TRINCULO",
        722: "FRANCISCO",
        723: "MARGARET",
        724: "FERDINAND",
        725: "PERDITA",
        726: "MAB",
        727: "CUPID",
        899: "NEPTUNE",
        801: "TRITON",
        802: "NEREID",
        803: "NAIAD",
        804: "THALASSA",
        805: "DESPINA",
        806: "GALATEA",
        807: "LARISSA",
        808: "PROTEUS",
        809: "HALIMEDE",
        810: "PSAMATHE",
        811: "SAO",
        812: "LAOMEDEIA",
        813: "NESO",
        999: "PLUTO",
        901: "CHARON",
        902: "NIX",
        903: "HYDRA",
        904: "KERBEROS",
        905: "STYX",
        -1: "GEOTAIL",
        -3: "MOM",
        -3: "MARS ORBITER MISSION",
        -5: "AKATSUKI",
        -5: "VCO",
        -5: "PLC",
        -5: "PLANET-C",
        -6: "P6",
        -6: "PIONEER-6",
        -7: "P7",
        -7: "PIONEER-7",
        -8: "WIND",
        -12: "VENUS ORBITER",
        -12: "P12",
        -12: "PIONEER 12",
        -12: "LADEE",
        -13: "POLAR",
        -18: "MGN",
        -18: "MAGELLAN",
        -18: "LCROSS",
        -20: "P8",
        -20: "PIONEER-8",
        -21: "SOHO",
        -23: "P10",
        -23: "PIONEER-10",
        -24: "P11",
        -24: "PIONEER-11",
        -25: "LP",
        -25: "LUNAR PROSPECTOR",
        -27: "VK1",
        -27: "VIKING 1 ORBITER",
        -28: "JUPITER ICY MOONS EXPLORER",
        -28: "JUICE",
        -29: "STARDUST",
        -29: "SDU",
        -29: "NEXT",
        -30: "VK2",
        -30: "VIKING 2 ORBITER",
        -30: "DS-1",
        -31: "VG1",
        -31: "VOYAGER 1",
        -32: "VG2",
        -32: "VOYAGER 2",
        -33: "NEOS",
        -33: "NEO SURVEYOR",
        -37: "HYB2",
        -37: "HAYABUSA 2",
        -37: "HAYABUSA2",
        -39: "LUNAR POLAR HYDROGEN MAPPER",
        -39: "LUNAH-MAP",
        -40: "CLEMENTINE",
        -41: "MEX",
        -41: "MARS EXPRESS",
        -43: "IMAP",
        -44: "BEAGLE2",
        -44: "BEAGLE 2",
        -45: "JNSA",
        -45: "JANUS_A",
        -46: "MS-T5",
        -46: "SAKIGAKE",
        -47: "PLANET-A",
        -47: "SUISEI",
        -47: "GNS",
        -47: "GENESIS",
        -48: "HUBBLE SPACE TELESCOPE",
        -48: "HST",
        -49: "LUCY",
        -53: "MARS PATHFINDER",
        -53: "MPF",
        -53: "MARS ODYSSEY",
        -53: "MARS SURVEYOR 01 ORBITER",
        -55: "ULYSSES",
        -57: "LUNAR ICECUBE",
        -58: "VSOP",
        -58: "HALCA",
        -59: "RADIOASTRON",
        -61: "JUNO",
        -62: "EMM",
        -62: "EMIRATES MARS MISSION",
        -64: "ORX",
        -64: "OSIRIS-REX",
        -65: "MCOA",
        -65: "MARCO-A",
        -66: "VEGA 1",
        -66: "MCOB",
        -66: "MARCO-B",
        -67: "VEGA 2",
        -68: "MERCURY MAGNETOSPHERIC ORBITER",
        -68: "MMO",
        -68: "BEPICOLOMBO MMO",
        -70: "DEEP IMPACT IMPACTOR SPACECRAFT",
        -72: "JNSB",
        -72: "JANUS_B",
        -74: "MRO",
        -74: "MARS RECON ORBITER",
        -76: "CURIOSITY",
        -76: "MSL",
        -76: "MARS SCIENCE LABORATORY",
        -77: "GLL",
        -77: "GALILEO ORBITER",
        -78: "GIOTTO",
        -79: "SPITZER",
        -79: "SPACE INFRARED TELESCOPE FACILITY",
        -79: "SIRTF",
        -81: "CASSINI ITL",
        -82: "CAS",
        -82: "CASSINI",
        -84: "PHOENIX",
        -85: "LRO",
        -85: "LUNAR RECON ORBITER",
        -85: "LUNAR RECONNAISSANCE ORBITER",
        -86: "CH1",
        -86: "CHANDRAYAAN-1",
        -90: "CASSINI SIMULATION",
        -93: "NEAR EARTH ASTEROID RENDEZVOUS",
        -93: "NEAR",
        -94: "MO",
        -94: "MARS OBSERVER",
        -94: "MGS",
        -94: "MARS GLOBAL SURVEYOR",
        -95: "MGS SIMULATION",
        -96: "PARKER SOLAR PROBE",
        -96: "SPP",
        -96: "SOLAR PROBE PLUS",
        -97: "TOPEX/POSEIDON",
        -98: "NEW HORIZONS",
        -107: "TROPICAL RAINFALL MEASURING MISSION",
        -107: "TRMM",
        -112: "ICE",
        -116: "MARS POLAR LANDER",
        -116: "MPL",
        -117: "EDL DEMONSTRATOR MODULE",
        -117: "EDM",
        -117: "EXOMARS 2016 EDM",
        -119: "MARS_ORBITER_MISSION_2",
        -119: "MOM2",
        -121: "MERCURY PLANETARY ORBITER",
        -121: "MPO",
        -121: "BEPICOLOMBO MPO",
        -127: "MARS CLIMATE ORBITER",
        -127: "MCO",
        -130: "MUSES-C",
        -130: "HAYABUSA",
        -131: "SELENE",
        -131: "KAGUYA",
        -135: "DART",
        -135: "DOUBLE ASTEROID REDIRECTION TEST",
        -140: "EPOCH",
        -140: "DIXI",
        -140: "EPOXI",
        -140: "DEEP IMPACT FLYBY SPACECRAFT",
        -142: "TERRA",
        -142: "EOS-AM1",
        -143: "TRACE GAS ORBITER",
        -143: "TGO",
        -143: "EXOMARS 2016 TGO",
        -144: "SOLO",
        -144: "SOLAR ORBITER",
        -146: "LUNAR-A",
        -148: "DFLY",
        -148: "DRAGONFLY",
        -150: "CASSINI PROBE",
        -150: "HUYGENS PROBE",
        -150: "CASP",
        -151: "AXAF",
        -151: "CHANDRA",
        -152: "CH2O",
        -152: "CHANDRAYAAN-2 ORBITER",
        -153: "CH2L",
        -153: "CHANDRAYAAN-2 LANDER",
        -154: "AQUA",
        -155: "KPLO",
        -155: "KOREAN PATHFINDER LUNAR ORBITER",
        -156: "ADITYA",
        -156: "ADIT",
        -159: "EURC",
        -159: "EUROPA CLIPPER",
        -164: "LUNAR FLASHLIGHT",
        -165: "MAP",
        -166: "IMAGE",
        -168: "PERSEVERANCE",
        -168: "MARS 2020",
        -168: "MARS2020",
        -168: "M2020",
        -170: "JWST",
        -170: "JAMES WEBB SPACE TELESCOPE",
        -172: "EXOMARS SCC",
        -173: "EXOMARS SP",
        -174: "EXOMARS ROVER",
        -177: "GRAIL-A",
        -178: "PLANET-B",
        -178: "NOZOMI",
        -181: "GRAIL-B",
        -183: "CLUSTER 1",
        -185: "CLUSTER 2",
        -188: "MUSES-B",
        -189: "NSYT",
        -189: "INSIGHT",
        -190: "SIM",
        -194: "CLUSTER 3",
        -196: "CLUSTER 4",
        -197: "EXOMARS_LARA",
        -197: "LARA",
        -198: "INTEGRAL",
        -198: "NASA-ISRO SAR MISSION",
        -198: "NISAR",
        -200: "CONTOUR",
        -202: "MAVEN",
        -203: "DAWN",
        -205: "SOIL MOISTURE ACTIVE AND PASSIVE",
        -205: "SMAP",
        -210: "LICIACUBE",
        -212: "STV51",
        -213: "STV52",
        -214: "STV53",
        -226: "ROSETTA",
        -227: "KEPLER",
        -228: "GLL PROBE",
        -228: "GALILEO PROBE",
        -234: "STEREO AHEAD",
        -235: "STEREO BEHIND",
        -236: "MESSENGER",
        -238: "SMART1",
        -238: "SM1",
        -238: "S1",
        -238: "SMART-1",
        -239: "MARTIAN MOONS EXPLORATION",
        -239: "MMX",
        -240: "SMART LANDER FOR INVESTIGATING MOON",
        -240: "SLIM",
        -242: "LUNAR TRAILBLAZER",
        -243: "VIPER",
        -248: "VEX",
        -248: "VENUS EXPRESS",
        -253: "OPPORTUNITY",
        -253: "MER-1",
        -254: "SPIRIT",
        -254: "MER-2",
        -255: "PSYC",
        -301: "HELIOS 1",
        -302: "HELIOS 2",
        -362: "RADIATION BELT STORM PROBE A",
        -363: "RADIATION BELT STORM PROBE B",
        -500: "SELENE Rstar",
        -502: "SELENE Vstar",
        -550: "MARS96",
        -652: "BEPICOLOMBO MTM",
        -750: "SPRINT-A",
        2000001: "CERES",
        2000002: "PALLAS",
        2000004: "VESTA",
        2000016: "PSYCHE",
        2000021: "LUTETIA",
        2000052: "52_EUROPA",
        2000052: "52 EUROPA",
        2000216: "KLEOPATRA",
        2000253: "MATHILDE",
        2000433: "EROS",
        2000511: "DAVIDA",
        2002867: "STEINS",
        2004015: "WILSON-HARRINGTON",
        2004179: "TOUTATIS",
        2009969: "1992KD",
        2009969: "BRAILLE",
        2025143: "ITOKAWA",
        2101955: "BENNU",
        2162173: "RYUGU",
        2431010: "IDA",
        2431011: "DACTYL",
        2486958: "ARROKOTH",
        9511010: "GASPRA",
        20000617: "PATROCLUS BARYCENTER",
        20003548: "EURYBATES BARYCENTER",
        20011351: "LEUCUS",
        20015094: "POLYMELE",
        20021900: "ORUS",
        20052246: "DONALDJOHANSON",
        20065803: "DIDYMOS BARYCENTER",
        120000617: "MENOETIUS",
        120003548: "QUETA",
        120065803: "DIMORPHOS",
        920000617: "PATROCLUS",
        920003548: "EURYBATES",
        920065803: "DIDYMOS",
    }
)

def naifid2name(naifid):
    """Retreives the planet/comet/spacecraft name from its NAIF id.

    Args:
        name (:class:`int`): The NAIF id.

    Returns:
        :class:`string` the name of the the planet/comet/spacecraft
    """
    if naifid in naif_dict:
        return naif_dict[naifid]
    else:
        return "Not Found"