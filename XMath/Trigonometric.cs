using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost
{
    public static partial class XMath
    {

        public static double sec(double x)
        {
            return 1/Math.Cos(x);
        }

        public static double cosec(double x)
        {
            return 1/Math.Sin(x);
        }

        public static double acot(double x)
        {
            return Math.PI/2 - Math.Atan(x);
        }

        public static double asec(double x)
        {
            if (x == 0) throw new ArgumentException("Asec is not defined for x = 0");
            return Math.Acos(1/x);
        }

        public static double acosec(double x)
        {
            if (x == 0) throw new ArgumentException("Acsc is not defined for x = 0");
            return Math.Asin(1 / x);
        }

    }
}
