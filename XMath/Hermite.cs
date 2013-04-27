using System;
using System.Collections.Generic;
using System.Text;

namespace CSBoost
{
    partial class XMath
    {

        public static double hermite(uint n, double x)
        {
            return hermite_imp(n, x);
        }

        public static double hermite_next(uint n, double x, double Hn, double Hnm1)
        {
           return (2 * x * Hn - 2 * n * Hnm1);
        }

        private static double hermite_imp(uint n, double x)
        {
           double p0 = 1;
           double p1 = 2 * x;

           if(n == 0) return p0;

           uint c = 1;

           while(c < n)
           {
              swap(ref p0, ref p1);
              p1 = hermite_next(c, x, p0, p1);
              ++c;
           }
           return p1;
        }
    }
}
