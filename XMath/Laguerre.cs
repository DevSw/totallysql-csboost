using System;
using System.Collections.Generic;
using System.Text;

namespace CSBoost
{
    partial class XMath
    {

        public static double laguerre(uint n, double x)
        {
            return laguerre_imp(n, x);
        }

        public static double laguerre(uint n, uint m, double x)
        {
            return laguerre_imp(n, m, x);
        }

        public static double laguerre_next(uint n, double x, double Ln, double Lnm1)
        {
           return ((2 * n + 1 - x) * Ln - n * Lnm1) / (n + 1);
        }

        public static double laguerre_next(uint n, uint l, double x, double Pl, double Plm1)
        {
            return ((2 * n + l + 1 - x) * Pl - (n + l) * Plm1) / (n + 1);
        }

        #region implementation

        private static double laguerre_imp(uint n, double x)
        {
           double p0 = 1;
           double p1 = 1 - x;

           if(n == 0)
              return p0;

           uint c = 1;

           while(c < n)
           {
              swap(ref p0, ref p1);
              p1 = laguerre_next(c, x, p0, p1);
              ++c;
           }
           return p1;
        }

        private static double laguerre_imp(uint n, uint m, double x)
        {
           // Special cases:
           if(m == 0) return laguerre(n, x);

           double p0 = 1;
           
           if(n == 0) return p0;

           double p1 = m + 1 - x;

           uint c = 1;

           while(c < n)
           {
              swap(ref p0, ref p1);
              p1 = laguerre_next(c, m, x, p0, p1);
              ++c;
           }
           return p1;
        }

        #endregion

    }
}
