using System;
using System.Collections.Generic;
using System.Text;

namespace CSBoost
{
    public partial class XMath
    {
        public static long Bin2Dec(string bin)
        {
            long result = Convert.ToInt64(bin, 2);
            return result;
        }
        public static long Hex2Dec(string hex)
        {
            long result = Convert.ToInt64(hex, 16);
            return result;
        }
        public static long Oct2Dec(string oct)
        {
            long result = Convert.ToInt64(oct, 8);
            return result;
        }
        public static long Base642Dec(string Base64)
        {
            byte[] ba = Convert.FromBase64String(Base64);
            //if (ba.Length > 4) throw new ArgumentException("Base642Dec: Input number is too large to be converted to a 64-bit integer");
            if (BitConverter.IsLittleEndian) Array.Reverse(ba);
            long result = BitConverter.ToInt64(ba, 0);
            return result;
        }
        public static string Dec2Bin(long dec)
        {
            return Convert.ToString(dec, 2);
        }
        public static string Dec2Hex(long dec)
        {
            return Convert.ToString(dec, 16);
        }
        public static string Dec2Oct(long dec)
        {
            return Convert.ToString(dec, 8);
        }
        public static string Dec2Base64(long dec)
        {
            byte[] ba = BitConverter.GetBytes(dec);
            if (BitConverter.IsLittleEndian) Array.Reverse(ba);
            return Convert.ToBase64String(ba);
        }
        public static string Oct2Bin(string oct)
        {
            return Dec2Bin(Oct2Dec(oct));
        }
        public static string Oct2Hex(string oct)
        {
            return Dec2Hex(Oct2Dec(oct));
        }
        public static string Oct2Base64(string oct)
        {
            return Dec2Base64(Oct2Dec(oct));
        }
        public static string Hex2Bin(string hex)
        {
            return Dec2Bin(Hex2Dec(hex));
        }
        public static string Hex2Oct(string hex)
        {
            return Dec2Oct(Hex2Dec(hex));
        }
        public static string Hex2Base64(string hex)
        {
            return Dec2Base64(Hex2Dec(hex));
        }
        public static string Bin2Hex(string bin)
        {
            return Dec2Hex(Bin2Dec(bin));
        }
        public static string Bin2Oct(string bin)
        {
            return Dec2Oct(Bin2Dec(bin));
        }
        public static string Bin2Base64(string bin)
        {
            return Dec2Base64(Bin2Dec(bin));
        }
        public static string Base642Hex(string Base64)
        {
            return Dec2Hex(Base642Dec(Base64));
        }
        public static string Base642Oct(string Base64)
        {
            return Dec2Oct(Base642Dec(Base64));
        }
        public static string Base642Bin(string Base64)
        {
            return Dec2Bin(Base642Dec(Base64));
        }

        #region conversion factors
        public class Unit
        {
            public string type, name, units, master;
            public double factor;
            public bool prefixable;

            public Unit(string t, string n, string u, double f, string m, bool p)
            {
                type = t;
                name = n;
                units = u;
                master = m;
                factor = f;
                prefixable = p;
            }

            public string ToString()
            {
                string result = string.Format("{0:G}: {1:G} [{2:G}] (= {3:G} {4:G})", type, units, name, factor, master);
                return result;
            }
        }

        public static readonly List<Unit> units = new List<Unit>(new Unit[]
            {
                new Unit("acceleration", "feet per hour per second", "fph/s", 0.00008466667, "m/s2", false),
                new Unit("acceleration", "feet per minute per second", "fpm/s", 0.00508, "m/s2", false),
                new Unit("acceleration", "feet per second squared", "fps2", 0.3048, "m/s2", false),
                new Unit("acceleration", "standard gravity", "G", 9.80665, "m/s2", false),
                new Unit("acceleration", "inches per minute per second", "ipm/s", 0.0004233333, "m/s2", false),
                new Unit("acceleration", "inches per second squared", "ips2", 0.0254, "m/s2", false),
                new Unit("acceleration", "knot per second", "kn/s", 0.5144444, "m/s2", false),
                new Unit("acceleration", "metres per second squared", "m/s2", 1, "m/s2", true),
                new Unit("acceleration", "miles per hour per second", "mph/s", 0.44704, "m/s2", false),
                new Unit("acceleration", "miles per minute per second", "mpm/s", 26.8224, "m/s2", false),
                new Unit("acceleration", "miles per second squared", "mps2", 1609.344, "m/s2", false),
                new Unit("angle", "minute (arc)", "arc min", 0.000290888, "rad", false),
                new Unit("angle", "second (arc)", "arc sec", 0.000004848137, "rad", false),
                new Unit("angle", "degree (arc)", "arc deg", 0.017453293, "rad", false),
                new Unit("angle", "gradian", "grad", 0.015707963, "rad", true),
                new Unit("angle", "radian", "rad", 1, "rad", false),
                new Unit("area", "acre", "ac", 4046.8564224, "m2", false),
                new Unit("area", "circular inch", "circ in", 0.0005067075, "m2", false),
                new Unit("area", "circular thou", "circ thou", 0.0000000005067075, "m2", false),
                new Unit("area", "hectare", "ha", 10000, "m2", false),
                new Unit("area", "metres squared", "m2", 1, "m2", true),
                new Unit("area", "square chain", "sq ch", 404.68564224, "m2", false),
                new Unit("area", "square foot", "sq ft", 0.09290304, "m2", false),
                new Unit("area", "square inch", "sq in", 0.00064516, "m2", false),
                new Unit("area", "square mile", "sq mi", 2589988.110336, "m2", false),
                new Unit("area", "square thou", "sq thou", 0.00000000064516, "m2", false),
                new Unit("area", "square yard", "sq yd", 0.83612736, "m2", false),
                new Unit("density", "grams per litre", "g/L", 1000, "g/m3", true),
                new Unit("density", "grams per cubic metre", "g/m3", 1, "g/m3", true),
                new Unit("density", "grams per millilitre", "g/mL", 1000000, "g/m3", true),
                new Unit("density", "pounds per cubic foot", "lb/ft3", 16018.46337, "g/m3", false),
                new Unit("density", "pounds per gallon (Imperial)", "lb/gal", 99776.37266, "g/m3", false),
                new Unit("density", "pounds per gallon (US)", "lb/gal", 119826.4273, "g/m3", false),
                new Unit("density", "pounds per cubic inch", "lb/in3", 27679904.71, "g/m3", false),
                new Unit("density", "ounces per cubic foot", "oz/ft3", 1001.153961, "g/m3", false),
                new Unit("density", "ounces per gallon (Imperial)", "oz/gal (Imp)", 6236.023291, "g/m3", false),
                new Unit("density", "ounces per gallon (US)", "oz/gal (US)", 7489.151707, "g/m3", false),
                new Unit("density", "ounces per cubic inch", "oz/in3", 1729994.044, "g/m3", false),
                new Unit("energy", "barrel of oil equivalent", "bboe", 6120000000, "J", false),
                new Unit("energy", "British thermal unit", "BTU", 1055.05585262, "J", false),
                new Unit("energy", "calorie", "cal", 4.1868, "J", false),
                new Unit("energy", "erg", "e", 0.0000001, "J", false),
                new Unit("energy", "electronvolt", "eV", 1.60217653E-19, "J", true),
                new Unit("energy", "foot-pound", "flb", 1.3558179483314, "J", false),
                new Unit("energy", "foot-poundal", "ft pdl", 0.0421401100938048, "J", false),
                new Unit("energy", "horsepower-hour", "HPh", 2684519.53769617, "J", false),
                new Unit("energy", "inch-pound force", "in lbf", 0.112984829027616, "J", false),
                new Unit("energy", "joule", "J", 1, "J", true),
                new Unit("energy", "kilocalorie", "kcal", 4186.8, "J", false),
                new Unit("energy", "ton of coal equivalent", "TCE", 29307600000, "J", false),
                new Unit("energy", "therm (E.C.)", "therm(EC)", 105505585.262, "J", false),
                new Unit("energy", "therm (U.S.)", "therm(US)", 105480400, "J", false),
                new Unit("energy", "ton of oil equivalent", "TOE", 41868000000, "J", false),
                new Unit("energy", "ton of TNT", "tTNT", 4184000000, "J", false),
                new Unit("energy", "watt-hour", "Wh", 3600, "J", true),
                new Unit("flow", "cubic feet per minute", "CFM", 0.0004719474432, "m3/s", false),
                new Unit("flow", "cubic feet per second", "ft3/s", 0.028316846592, "m3/s", false),
                new Unit("flow", "gallons (US) per day", "GPD", 0.00000004381263638, "m3/s", false),
                new Unit("flow", "gallons (US) per hour", "GPH", 0.000001051503273, "m3/s", false),
                new Unit("flow", "gallons (US) per minute", "GPM", 0.0000630901964, "m3/s", false),
                new Unit("flow", "cubic inches per minute", "in3/min", 0.00000027311773, "m3/s", false),
                new Unit("flow", "cubic inches per second", "in3/s", 0.000016387064, "m3/s", false),
                new Unit("flow", "litres per minute", "l/min", 0.000016, "m3/s", true),
                new Unit("flow", "litres per second", "l/s", 0.000016, "m3/s", true),
                new Unit("flow", "metres cubed per second", "m3/s", 1, "m3/s", true),
                new Unit("force", "dyne", "dyn", 0.00001, "N", false),
                new Unit("force", "kilogram-force", "kgf", 9.80665, "N", false),
                new Unit("force", "kip", "kip", 4448.2216152605, "N", false),
                new Unit("force", "pound-force", "lbf", 4.4482216152605, "N", false),
                new Unit("force", "newtons", "N", 1, "N", false),
                new Unit("force", "ounce-force", "ozf", 0.278013850953781, "N", false),
                new Unit("force", "poundal", "pdl", 0.138254954376, "N", false),
                new Unit("force", "ton-force", "tnf", 8896.443230521, "N", false),
                new Unit("length", "micron", "u", 0.000001, "m", false),
                new Unit("length", "ångström", "ang", 0.0000000001, "m", false),
                new Unit("length", "astronomical unit", "AU", 149597871464, "m", false),
                new Unit("length", "cable length", "cable", 185.2, "m", false),
                new Unit("length", "chain", "chain", 20.11684, "m", false),
                new Unit("length", "fathom", "fm", 1.8288, "m", false),
                new Unit("length", "foot", "ft", 0.3048, "m", false),
                new Unit("length", "furlong", "fur", 201.168, "m", false),
                new Unit("length", "hand", "hand", 0.1016, "m", false),
                new Unit("length", "inch", "in", 0.0254, "m", false),
                new Unit("length", "league", "lea", 4828.032, "m", false),
                new Unit("length", "light-year", "ly", 9460730472580800, "m", false),
                new Unit("length", "metre", "m", 1, "m", true),
                new Unit("length", "mile", "mi", 1609.344, "m", false),
                new Unit("length", "nautical mile", "Nmi", 1852, "m", false),
                new Unit("length", "pace", "pace", 0.762, "m", false),
                new Unit("length", "parsec", "pc", 30856778200000000, "m", false),
                new Unit("length", "pica", "pica", 0.004233333333, "m", false),
                new Unit("length", "point", "pt",  0.000352777778, "m", false),
                new Unit("length", "thou", "thou", 0.0000254, "m", false),
                new Unit("length", "twip", "twp", 0.000017638, "m", false),
                new Unit("length", "Yard", "yd", 0.9144, "m", false),
                new Unit("magnetic field", "gauss", "ga", 0.0001, "T", false),
                new Unit("magnetic field", "tesla", "T", 1, "T", true),
                new Unit("mass", "atomic mass unit", "u", 1.66053872E-24, "g", false),
                new Unit("mass", "carat (metric)", "ct", 0.2, "g", false),
                new Unit("mass", "hundredweight (long)", "cwt", 50802.34544, "g", false),
                new Unit("mass", "dram", "dr", 1.7718451953125, "g", false),
                new Unit("mass", "pennyweight", "dwt", 1.55517384, "g", false),
                new Unit("mass", "electronvolt", "eVm", 1.7826E-33, "g", false),
                new Unit("mass", "gram", "g", 1, "g", true),
                new Unit("mass", "grain", "gr", 0.06479891, "g", false),
                new Unit("mass", "carat", "kt", 0.205196548333, "g", false),
                new Unit("mass", "pound (avoirdupois)", "lb av", 453.59237, "g", false),
                new Unit("mass", "pound (avoirdupois)", "lbm", 453.59237, "g", false),
                new Unit("mass", "pound (troy)", "lb t", 373.2417216, "g", false),
                new Unit("mass", "electron rest mass", "me", 9.10938215E-28, "g", false),
                new Unit("mass", "ounce (avoirdupois)", "oz av", 28.349523125, "g", false),
                new Unit("mass", "ounce (avoirdupois)", "ozm", 28.349523125, "g", false),
                new Unit("mass", "ounce (troy)", "oz t", 31.1034768, "g", false),
                new Unit("mass", "hundredweight (short)", "sh cwt", 45359.237, "g", false),
                new Unit("mass", "ton, short", "sh tn", 907184.74, "g", false),
                new Unit("mass", "slug", "sg", 14593.903, "g", false),
                new Unit("mass", "stone", "st", 6350.29318, "g", false),
                new Unit("mass", "tonne", "t", 1000000, "g", true),
                new Unit("mass", "ton, long", "ton", 1016046.9088, "g", false),
                new Unit("power", "BTU per hour", "BTU/h", 0.293071, "W", false),
                new Unit("power", "BTU per minute", "BTU/min", 17.584264, "W", false),
                new Unit("power", "BTU per second", "BTU/s", 1055.05585262, "W", false),
                new Unit("power", "calorie per second", "cal/s", 4.1868, "W", false),
                new Unit("power", "foot-pound-force per hour", "ft lbf/h", 0.0003766161, "W", false),
                new Unit("power", "foot-pound-force per minute", "ft lbf/min", 0.0225969658055233, "W", false),
                new Unit("power", "foot-pound-force per second", "ft lbf/s", 1.3558179483314, "W", false),
                new Unit("power", "horsepower", "HP", 745.69987158227, "W", false),
                new Unit("power", "watt", "W", 1, "W", true),
                new Unit("pressure", "atmosphere (standard)", "atm", 101325, "Pa", false),
                new Unit("pressure", "bar", "bar", 100000, "Pa", false),
                new Unit("pressure", "metre of water", "mH2O", 9806.38, "Pa", true),
                new Unit("pressure", "metre of mercury", "mHg", 133322.368421, "Pa", true),
                new Unit("pressure", "foot of water", "ftH2O", 2988.98, "Pa", false),
                new Unit("pressure", "foot of mercury", "ftHg", 40636.66, "Pa", false),
                new Unit("pressure", "inch of water", "inH2O", 249.082, "Pa", false),
                new Unit("pressure", "inch of mercury", "inHg", 3386.389, "Pa", false),
                new Unit("pressure", "gram-force per square millimetre", "gf/mm2", 9806.65, "Pa", true),
                new Unit("pressure", "kip per square inch", "ksi", 6894757, "Pa", false),
                new Unit("pressure", "pascal", "Pa", 1, "Pa", true),
                new Unit("pressure", "pascal", "p", 1, "Pa", true),
                new Unit("pressure", "poundal per square foot", "pdl/sq ft", 1.488164, "Pa", false),
                new Unit("pressure", "pound per square foot", "psf", 47.88026, "Pa", false),
                new Unit("pressure", "pound per square inch", "psi", 6894.757, "Pa", false),
                new Unit("pressure", "torr", "torr", 133.3223494, "Pa", false),
                new Unit("pressure", "short ton per square foot", "sh tn/sq ft", 95760.518, "Pa", false),
                new Unit("speed", "speed of light in vacuum", "c", 299792458, "m/s", false),
                new Unit("speed", "feet per hour", "fph", 0.00008466667, "m/s", false),
                new Unit("speed", "feet per minute", "fpm", 0.00508, "m/s", false),
                new Unit("speed", "feet per second", "fps", 0.3048, "m/s", false),
                new Unit("speed", "inches per minute", "ipm", 0.000423333, "m/s", false),
                new Unit("speed", "inches per second", "ips", 0.0254, "m/s", false),
                new Unit("speed", "knot", "kn", 0.514444, "m/s", false),
                new Unit("speed", "metres per hour", "m/h", 0.0002777778, "m/s", true),
                new Unit("speed", "metres per second", "m/s", 1, "m/s", true),
                new Unit("speed", "miles per hour", "mph", 0.44704, "m/s", false),
                new Unit("speed", "miles per minute", "mpm", 26.8224, "m/s", false),
                new Unit("speed", "miles per second", "mps", 1609.344, "m/s", false),
                new Unit("temperature", "degree Celsius", "C", 1, "K", false),
                new Unit("temperature", "degree Fahrenheit", "F", 0.555555555555556, "K", false),
                new Unit("temperature", "kelvin", "K", 1, "K", false),
                new Unit("time", "Day", "day", 86400, "sec", false),
                new Unit("time", "Hour", "hr", 3600, "sec", false),
                new Unit("time", "Minute", "mn", 60, "sec", false),
                new Unit("time", "Month", "mo", 2629700, "sec", false),
                new Unit("time", "Second", "sec", 1, "sec", true),
                new Unit("time", "Week", "wk", 604800, "sec", false),
                new Unit("time", "Year", "yr", 31557600, "sec", false),
                new Unit("volume", "acre-foot", "ac ft", 1233.48183754752, "m3", false),
                new Unit("volume", "barrel (petroleum)", "bbl", 0.158987294928, "m3", false),
                new Unit("volume", "barrel (Imperial)", "bl (Imp)", 0.16365924, "m3", false),
                new Unit("volume", "barrel (US dry)", "bl (US)", 0.115628198985075, "m3", false),
                new Unit("volume", "cubic foot", "cu ft", 0.028316846592, "m3", false),
                new Unit("volume", "cubic inch", "cu in", 0.000016387064, "m3", false),
                new Unit("volume", "cubic mile", "cu mi", 4168181825.44057, "m3", false),
                new Unit("volume", "cubic yard", "cu yd", 0.764554857984, "m3", false),
                new Unit("volume", "cup (metric)", "cup", 0.00025, "m3", false),
                new Unit("volume", "cup (US)", "cup (US)", 0.00024, "m3", false),
                new Unit("volume", "ton (displacement)", "dt", 0.99108963072, "m3", false),
                new Unit("volume", "barrel (US fluid)", "fl bl (US)", 0.119240471196, "m3", false),
                new Unit("volume", "fluid ounce (Imperial)", "fl oz (Imp)", 0.0000284130625, "m3", false),
                new Unit("volume", "fluid ounce (US)", "fl oz (US)", 0.00003, "m3", false),
                new Unit("volume", "gallon (Imperial)", "gal (Imp)", 0.00454609, "m3", false),
                new Unit("volume", "gallon (US)", "gal (US)", 0.003785411784, "m3", false),
                new Unit("volume", "gill (Imperial)", "gi (Imp)", 0.0001420653125, "m3", false),
                new Unit("volume", "gill (US)", "gi (US)", 0.00011829411825, "m3", false),
                new Unit("volume", "ton (register)", "GRT", 2.8316846592, "m3", false),
                new Unit("volume", "litre", "l", 0.001, "m3", true),
                new Unit("volume", "metre cubed", "m3", 1, "m3", true),
                new Unit("volume", "cubic centimetre", "cc", 0.000001, "m3", false),
                new Unit("volume", "ton (freight)", "MT", 1.13267386368, "m3", false),
                new Unit("volume", "pint (Imperial)", "pt (Imp)", 0.00056826125, "m3", false),
                new Unit("volume", "pint (US)", "pt (US)", 0.000473176473, "m3", false),
                new Unit("volume", "quart (Imperial)", "qt (Imp)", 0.0011365225, "m3", false),
                new Unit("volume", "quart (US)", "qt (US)", 0.000946352946, "m3", false),
                new Unit("volume", "tablespoon", "tbs", 0.000015, "m3", false),
                new Unit("volume", "teaspoon", "tsp", 0.000005, "m3", false)
            });

        class UnitPrefix
        {
            public string name, prefix;
            public double factor;
            public UnitPrefix(string n, string p, double f)
            {
                name = n;
                prefix = p;
                factor = f;
            }
        }

        static readonly List<UnitPrefix> unitprefixes = new List<UnitPrefix>(new UnitPrefix[]
            {
                new UnitPrefix("exa", "E", 1e18),
                new UnitPrefix("peta", "P", 1e15),
                new UnitPrefix("tera", "T", 1e12),
                new UnitPrefix("giga", "G", 1e9),
                new UnitPrefix("mega", "M", 1e6),
                new UnitPrefix("kilo", "k", 1e3),
                new UnitPrefix("hecto", "h", 1e2),
                new UnitPrefix("dekao", "e", 1e1),
                new UnitPrefix("deci", "d", 1e-1),
                new UnitPrefix("centi", "c", 1e-2),
                new UnitPrefix("milli", "m", 1e-3),
                new UnitPrefix("micro", "u", 1e-6),
                new UnitPrefix("nano", "n", 1e-9),
                new UnitPrefix("pico", "p", 1e-12),
                new UnitPrefix("femto", "f", 1e-15),
                new UnitPrefix("atto", "a", 1e-18)
            });

        #endregion

        public static double UConvert(double value, string from_unit, string to_unit)
        {
            if(from_unit == to_unit) return value;
            Unit FromU = null, ToU = null;
            UnitPrefix FromP = null, ToP = null;
            string type = "";
            foreach (Unit u in units)
            {
                if (u.units == from_unit) FromU = u;
                if (u.units == to_unit) ToU = u;
                if (FromU != null && ToU != null) break;
            }
            if(FromU == null || ToU == null)
            {
                if (FromU != null) type = FromU.type;
                else if (ToU != null) type = ToU.type;
                foreach(Unit u in units)
                {
                    if (u.prefixable && (type == u.type || type == ""))
                    {
                        foreach (UnitPrefix up in unitprefixes)
                        {
                            if (FromU == null && up.prefix + u.units == from_unit)
                            {
                                FromU = u;
                                FromP = up;
                            }

                            if (ToU == null && up.prefix + u.units == to_unit)
                            {
                                ToU = u;
                                ToP = up;
                            }
                        }
                    }
                    if (FromU != null && ToU != null) break;
                }
            }
            if (FromU == null) throw new ArgumentException(string.Format("Convert: 'From' units ({0:G}) not recognized - use SELECT * FROM [schema].ConversionUnits to show valid units", from_unit));
            if (ToU == null) throw new ArgumentException(string.Format("Convert: 'To' units ({0:G}) not recognized - use SELECT * FROM [schema].ConversionUnits to show valid units", to_unit));
            if (FromU.type != ToU.type) throw new ArgumentException(string.Format("Convert: 'To' units ({0:G} - of type {1:G}) and 'From' units ({2:G} - of type {3:G}) are not of the same type - conversion between them is not possible", to_unit, ToU.type, from_unit, FromU.type));
            double fromfactor = FromP == null ? FromU.factor : FromU.factor * FromP.factor;
            double tofactor = ToP == null ? ToU.factor : ToU.factor * ToP.factor;
            double result = value * fromfactor / tofactor;
            if (FromU.type == "temperature")
            {
                switch (FromU.units)
                {
                    case "K":
                        if (ToU.units == "C") result = value -273.15;
                        if (ToU.units == "F") result -= 459.67;
                        break;
                    case "C":
                        if (ToU.units == "K") result = value + 273.15;
                        if (ToU.units == "F") result += 32;
                        break;
                    case "F":
                        if (ToU.units == "C") result = (value - 32) * (5.0/9.0);
                        if (ToU.units == "K") result = (value + 459.67) * (5.0 / 9.0);
                        break;
                }
            }
            return result;
        }
    }
}
