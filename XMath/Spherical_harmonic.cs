using System;
using System.Collections.Generic;
using System.Text;

namespace CSBoost
{
    partial class XMath
    {
        public static double spherical_harmonic_r(uint n, int m, double theta, double phi)
        {
           bool sign = false;
           if(m < 0)
           {
              // Reflect and adjust sign if m < 0:
              sign = (m&1)>0;
              m = Math.Abs(m);
           }
           if((m&1)>0)
           {
              // Check phase if theta is outside [0, PI]:
              double mod = fmod(theta, 2 * Math.PI);
              if(mod < 0)
                 mod += 2 * Math.PI;
              if(mod > Math.PI)
                 sign = !sign;
           }
           // Get the value and adjust sign as required:
           double prefix = spherical_harmonic_prefix(n, m, theta);
           prefix *= Math.Cos(m * phi);
           return sign ? -prefix : prefix;
        }

        public static double spherical_harmonic_i(uint n, int m, double theta, double phi)
        {
           bool sign = false;
           if(m < 0)
           {
              // Reflect and adjust sign if m < 0:
              sign = !((m&1)>0);
              m = Math.Abs(m);
           }
           if((m&1)>0)
           {
              // Check phase if theta is outside [0, PI]:
              double mod = fmod(theta, 2 * Math.PI);
              if(mod < 0)
                 mod += 2 * Math.PI;
              if(mod > Math.PI)
                 sign = !sign;
           }
           // Get the value and adjust sign as required:
           double prefix = spherical_harmonic_prefix(n, m, theta);
           prefix *= Math.Sin(m * phi);
           return sign ? -prefix : prefix;
        }

        public static Complex spherical_harmonic(uint n, int m, double theta, double phi)
        {
           //
           // Sort out the signs:
           //
           bool r_sign = false;
           bool i_sign = false;
           if(m < 0)
           {
              // Reflect and adjust sign if m < 0:
              r_sign = (m&1)>0;
              i_sign = !r_sign;
              m = Math.Abs(m);
           }
           if((m&1)>0)
           {
              // Check phase if theta is outside [0, PI]:
              double mod = fmod(theta, 2 * Math.PI);
              if(mod < 0)
                 mod += 2 * Math.PI;
              if(mod > Math.PI)
              {
                 r_sign = !r_sign;
                 i_sign = !i_sign;
              }
           }
           //
           // Calculate the value:
           //
           double prefix = spherical_harmonic_prefix(n, m, theta);
           double r = prefix * Math.Cos(m * phi);
           double i = prefix * Math.Sin(m * phi);
           //
           // Add in the signs:
           //
           if(r_sign)
              r = -r;
           if(i_sign)
              i = -i;
           return new Complex(r, i);
        }

        private static double spherical_harmonic_prefix(uint n, int m, double theta)
        {
            if (m > n) return 0;

            double sin_theta = Math.Sin(theta);
            double x = Math.Cos(theta);

            double leg = legendre_p_imp((int)n, m, x, Math.Pow(Math.Abs(sin_theta), m));

            double prefix = gamma_delta_ratio(n - m + 1, 2 * m);
            prefix *= (2 * n + 1) / (4 * Math.PI);
            prefix = Math.Sqrt(prefix);
            return prefix * leg;
        }
    }
}
