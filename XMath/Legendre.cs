using System;
using System.Collections.Generic;
using System.Text;

namespace CSBoost
{
    partial class XMath
    {
        public static double legendre_next(uint l, double x, double Pl, double Plm1)
        {
           return ((2 * l + 1.0) * x * Pl - l * Plm1) / (l + 1.0);
        }

        public static double legendre_p(int l, double x)
        {
           if(l < 0) return legendre_imp((uint)-l-1, x, false);
           return legendre_imp((uint)l, x, false);
        }

        public static double legendre_q(uint l, double x)
        {
           return legendre_imp(l, x, true);
        }

        public static double legendre_next(uint l, uint m, double x, double Pl, double Plm1)
        {
           return ((2 * l + 1.0) * x * Pl - (l + m) * Plm1) / (l + 1.0 - m);
        }

        public static double legendre_p(int l, int m, double x)
        {
            return legendre_p_imp(l, m, x);
        }

        private static double legendre_p_imp(int l, int m, double x, double sin_theta_power)
        {
           // Error handling:
           if((x < -1) || (x > 1))
               throw new ArgumentException(string.Format("The associated Legendre Polynomial is defined for -1 <= x <= 1, but got x = {0:G}.", x));
           // Handle negative arguments first:
           if(l < 0) return legendre_p_imp(-l-1, m, x, sin_theta_power);
           if(m < 0)
           {
              int sign = (m&1) > 0 ? -1 : 1;
              return sign * gamma_ratio(l+m+1, l+1-m) * legendre_p_imp(l, -m, x, sin_theta_power);
           }
           // Special cases:
           if(m > l) return 0;
           if(m == 0) return legendre_p(l, x);

           double p0 = double_factorial((uint)(2 * m - 1)) * sin_theta_power;

           if((m&1) > 0) p0 *= -1;
           if(m == l) return p0;

           double p1 = x * (2 * m + 1) * p0;

           int n = m + 1;

           while(n < l)
           {
              swap(ref p0, ref p1);
              p1 = legendre_next((uint)n, (uint)m, x, p0, p1);
              ++n;
           }
           return p1;
        }

        private static double legendre_p_imp(int l, int m, double x)
        {
           return legendre_p_imp(l, m, x, Math.Pow(1 - x*x, Math.Abs(m)/2.0));
        }

        private static double legendre_imp(uint l, double x, bool second)
        {
           if((x < -1) || (x > 1)) throw new ArgumentException(string.Format("Legendre Polynomial: not defined for x < -1 or x > 1 (got {0:G}).", x));

           double p0, p1;
           if(second)
           {
              // A solution of the second kind (Q):
              p0 = (log1p(x) - log1p(-x)) / 2;
              p1 = x * p0 - 1;
           }
           else
           {
              // A solution of the first kind (P):
              p0 = 1;
              p1 = x;
           }
           if(l == 0)
              return p0;

           uint n = 1;

           while(n < l)
           {
              swap(ref p0, ref p1);
              p1 = legendre_next(n, x, p0, p1);
              ++n;
           }
           return p1;
        }
    }
}
