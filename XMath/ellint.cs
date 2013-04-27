using System;
using System.Collections.Generic;
using System.Text;

namespace CSBoost
{
    partial class XMath
    {
        public static double ellint_rc(double x, double y)
        {
            return ellint_rc_imp(x, y);
        }

        public static double ellint_rd(double x, double y, double z)
        {
            return ellint_rd_imp(x, y, z);
        }

        public static double ellint_rf(double x, double y, double z)
        {
            return ellint_rf_imp(x, y, z);
        }

        public static double ellint_rj(double x, double y, double z, double p)
        {
            return ellint_rj_imp(x, y, z, p);
        }

        public static double ellint_1(double k)
        {
           return ellint_k_imp(k);
        }

        public static double ellint_1(double k, double phi)
        {
           return ellint_f_imp(phi, k);
        }

        public static double ellint_2(double k)
        {
           return ellint_e_imp(k);
        }

        public static double ellint_2(double k, double phi)
        {
           return ellint_e_imp(phi, k);
        }

        public static double ellint_3(double k, double v)
        {
            return ellint_pi_imp(v, k, 1 - v);
        }

        public static double ellint_3(double k, double v, double phi)
        {
            return ellint_pi_imp(v, phi, k, 1 - v);
        }

        private static double ellint_rc_imp(double x, double y)
        {
            double value, S, u, lambda, tolerance, prefix;
            ulong k;

            if(x < 0) throw new Exception(string.Format("Argument x must be non-negative but got {0:G}", x));
            if(y == 0) throw new Exception(string.Format("Argument y must not be zero but got {0:G}", y));

            // error scales as the 6th power of tolerance
            tolerance = Math.Pow(4 * epsilon, 1.0 / 6);

            // for y < 0, the integral is singular, return Cauchy principal value
            if (y < 0)
            {
                prefix = Math.Sqrt(x / (x - y));
                x = x - y;
                y = -y;
            }
            else
               prefix = 1;

            // duplication:
            k = 1;
            do
            {
                u = (x + y + y) / 3;
                S = y / u - 1;               // 1 - x / u = 2 * S

                if (2 * Math.Abs(S) < tolerance) 
                   break;

                double sx = Math.Sqrt(x);
                double sy = Math.Sqrt(y);
                lambda = 2 * sx * sy + y;
                x = (x + lambda) / 4;
                y = (y + lambda) / 4;
                ++k;
            }while(k < max_series_iterations);
            check_series_iterations(k);

            // Taylor series expansion to the 5th order
            value = (1 + S * S * (3.0 / 10 + S * (1.0 / 7 + S * (3.0 / 8 + S * 9.0 / 22)))) / Math.Sqrt(u);

            return value * prefix;
        }

        private static double ellint_rd_imp(double x, double y, double z)
        {
            double value, u, lambda, sigma, factor, tolerance;
            double X, Y, Z, EA, EB, EC, ED, EE, S1, S2;
            ulong k;

            if (x < 0)
            {
               throw new Exception(string.Format("Argument x must be >= 0, but got {0:G}", x));
            }
            if (y < 0)
            {
               throw new Exception(string.Format("Argument y must be >= 0, but got {0:G}", y));
            }
            if (z <= 0)
            {
               throw new Exception(string.Format("Argument z must be > 0, but got {0:G}", z));
            }
            if (x + y == 0)
            {
               throw new Exception(string.Format("At most one argument can be zero, but got, x + y = {0:G}", x+y));
            }

            // error scales as the 6th power of tolerance
            tolerance = Math.Pow(epsilon / 3, 1.0/6);

            // duplication
            sigma = 0;
            factor = 1;
            k = 1;
            do
            {
                u = (x + y + z + z + z) / 5;
                X = (u - x) / u;
                Y = (u - y) / u;
                Z = (u - z) / u;
                if (Math.Max(Math.Abs(X), Math.Max(Math.Abs(Y), Math.Abs(Z))) < tolerance) 
                   break;
                double sx = Math.Sqrt(x);
                double sy = Math.Sqrt(y);
                double sz = Math.Sqrt(z);
                lambda = sy * (sx + sz) + sz * sx; //Math.Sqrt(x * y) + Math.Sqrt(y * z) + Math.Sqrt(z * x);
                sigma += factor / (sz * (z + lambda));
                factor /= 4;
                x = (x + lambda) / 4;
                y = (y + lambda) / 4;
                z = (z + lambda) / 4;
                ++k;
            }
            while(k < max_series_iterations);

            // Check to see if we gave up too soon:
            check_series_iterations(k);

            // Taylor series expansion to the 5th order
            EA = X * Y;
            EB = Z * Z;
            EC = EA - EB;
            ED = EA - 6 * EB;
            EE = ED + EC + EC;
            S1 = ED * (ED * 9.0 / 88 - Z * EE * 9.0 / 52 - 3.0 / 14);
            S2 = Z * (EE / 6 + Z * (-EC * 9.0 / 22 + Z * EA * 3.0 / 26));
            value = 3 * sigma + factor * (1 + S1 + S2) / (u * Math.Sqrt(u));

            return value;
        }

        private static double ellint_rf_imp(double x, double y, double z)
        {
            double value, X, Y, Z, E2, E3, u, lambda, tolerance;
            ulong k;

            if (x < 0 || y < 0 || z < 0)
            {
               throw new Exception(string.Format("Domain error, all arguments must be non-negative, only sensible result is {0:G}.", double.NaN));
            }
            if (x + y == 0 || y + z == 0 || z + x == 0)
            {
               throw new Exception(string.Format("Domain error, at most one argument can be zero, only sensible result is {0:G}.", double.NaN));
            }

            tolerance = Math.Pow(4*epsilon, 1.0/6);

            // duplication
            k = 1;
            do
            {
                u = (x + y + z) / 3;
                X = (u - x) / u;
                Y = (u - y) / u;
                Z = (u - z) / u;

                // Termination condition: 
                if (Math.Max(Math.Abs(X), Math.Max(Math.Abs(Y), Math.Abs(Z))) < tolerance) 
                   break; 

                double sx = Math.Sqrt(x);
                double sy = Math.Sqrt(y);
                double sz = Math.Sqrt(z);
                lambda = sy * (sx + sz) + sz * sx;
                x = (x + lambda) / 4;
                y = (y + lambda) / 4;
                z = (z + lambda) / 4;
                ++k;
            }
            while(k < max_series_iterations);

            // Check to see if we gave up too soon:
            check_series_iterations(k);

            // Taylor series expansion to the 5th order
            E2 = X * Y - Z * Z;
            E3 = X * Y * Z;
            value = (1 + E2*(E2/24 - E3*3.0/44 - 0.1) + E3/14) / Math.Sqrt(u);

            return value;
        }

        private static double ellint_rj_imp(double x, double y, double z, double p)
        {
            double value, u, lambda, alpha, beta, sigma, factor, tolerance;
            double X, Y, Z, P, EA, EB, EC, E2, E3, S1, S2, S3;
            ulong k;

            if (x < 0)
            {
               throw new Exception(string.Format("Argument x must be non-negative, but got x = {0:G}", x));
            }
            if(y < 0)
            {
               throw new Exception(string.Format("Argument y must be non-negative, but got y = {0:G}", y));
            }
            if(z < 0)
            {
               throw new Exception(string.Format("Argument z must be non-negative, but got z = {0:G}", z));
            }
            if(p == 0)
            {
               throw new Exception(string.Format("Argument p must not be zero, but got p = {0:G}", p));
            }
            if (x + y == 0 || y + z == 0 || z + x == 0)
            {
               throw new Exception(string.Format("At most one argument can be zero, only possible result is {0:G}.", double.NaN));
            }

            // error scales as the 6th power of tolerance
            tolerance = Math.Pow(1.0 * epsilon / 3, 1.0 / 6);

            // for p < 0, the integral is singular, return Cauchy principal value
            if (p < 0)
            {
               //
               // We must ensure that (z - y) * (y - x) is positive.
               // Since the integral is symmetrical in x, y and z
               // we can just permute the values:
               //
               if(x > y)
                  swap(ref x, ref y);
               if(y > z)
                  swap(ref y, ref z);
               if(x > y)
                  swap(ref x, ref y);

               double q = -p;
               double pmy = (z - y) * (y - x) / (y + q);  // p - y
                
               double p2 = pmy + y;
               value = ellint_rj(x, y, z, p2);
               value *= pmy;
               value -= 3 * ellint_rf(x, y, z);
               value += 3 * Math.Sqrt((x * y * z) / (x * z + p2 * q)) * ellint_rc(x * z + p2 * q, p2 * q);
               value /= (y + q);
               return value;
            }

            // duplication
            sigma = 0;
            factor = 1;
            k = 1;
            do
            {
                u = (x + y + z + p + p) / 5;
                X = (u - x) / u;
                Y = (u - y) / u;
                Z = (u - z) / u;
                P = (u - p) / u;
                
                if (Math.Max(Math.Max(Math.Abs(X), Math.Abs(Y)), Math.Max(Math.Abs(Z), Math.Abs(P))) < tolerance) 
                   break;

                double sx = Math.Sqrt(x);
                double sy = Math.Sqrt(y);
                double sz = Math.Sqrt(z);
                
                lambda = sy * (sx + sz) + sz * sx;
                alpha = p * (sx + sy + sz) + sx * sy * sz;
                alpha *= alpha;
                beta = p * (p + lambda) * (p + lambda);
                sigma += factor * ellint_rc(alpha, beta);
                factor /= 4;
                x = (x + lambda) / 4;
                y = (y + lambda) / 4;
                z = (z + lambda) / 4;
                p = (p + lambda) / 4;
                ++k;
            }
            while(k < max_series_iterations);

            // Check to see if we gave up too soon:
            check_series_iterations(k);

            // Taylor series expansion to the 5th order
            EA = X * Y + Y * Z + Z * X;
            EB = X * Y * Z;
            EC = P * P;
            E2 = EA - 3 * EC;
            E3 = EB + 2 * P * (EA - EC);
            S1 = 1 + E2 * (E2 * 9.0 / 88 - E3 * 9.0 / 52 - 3.0 / 14);
            S2 = EB * (1.0 / 6 + P * (-6.0 / 22 + P * 3.0 / 26));
            S3 = P * ((EA - EC) / 3 - P * EA * 3.0 / 22);
            value = 3 * sigma + factor * (S1 + S2 + S3) / (u * Math.Sqrt(u));

            return value;
        }

        private static double ellint_f_imp(double phi, double k)
        {

            if (Math.Abs(k) > 1)
            {
               throw new Exception(string.Format("Got k = {0:G}, function requires |k| <= 1", k));
            }

            bool invert = false;
            if(phi < 0)
            {
               phi = Math.Abs(phi);
               invert = true;
            }

            double result;

            if(phi >= double.MaxValue)
            {
                throw new OverflowException();
            }
            else if(phi > 1 / epsilon)
            {
               result = 2 * phi * ellint_k_imp(k / Math.PI);
            }
            else
            {
               double rphi = fmod(phi, (Math.PI / 2));
               double m = Math.Floor((2 * phi) / Math.PI);
               int s = 1;
               if(fmod(m, 2.0) > 0.5)
               {
                  m += 1;
                  s = -1;
                  rphi = Math.PI / 2 - rphi;
               }
               double sinp = Math.Sin(rphi);
               double cosp = Math.Cos(rphi);
               result = s * sinp * ellint_rf_imp((cosp * cosp), (1 - k * k * sinp * sinp), 1.0);
               if(m != 0)
               {
                  result += m * ellint_k_imp(k);
               }
            }
            return invert ? -result : result;
        }

        private static double ellint_k_imp(double k)
        {
            if (Math.Abs(k) > 1)
            {
               throw new Exception(string.Format("Got k = {0:G}, function requires |k| <= 1", k));
            }
            if (Math.Abs(k) == 1) throw new OverflowException();

            double x = 0;
            double y = 1 - k * k;
            double z = 1;
            double value = ellint_rf_imp(x, y, z);

            return value;
        }

        private static double ellint_e_imp(double phi, double k)
        {
            bool invert = false;
            if(phi < 0)
            {
               phi = Math.Abs(phi);
               invert = true;
            }

            double result;

            if(phi >= double.MaxValue) throw new OverflowException();
            else if(phi > 1 / epsilon)
            {
               result = 2 * phi * ellint_e_imp(k) / Math.PI;
            }
            else
            {
               double rphi = fmod(phi, Math.PI / 2);
               double m = Math.Floor((2 * phi) / Math.PI);
               int s = 1;
               if(fmod(m, 2.0) > 0.5)
               {
                  m += 1;
                  s = -1;
                  rphi = Math.PI / 2 - rphi;
               }
               double sinp = Math.Sin(rphi);
               double cosp = Math.Cos(rphi);
               double x = cosp * cosp;
               double t = k * k * sinp * sinp;
               double y = 1 - t;
               double z = 1;
               result = s * sinp * (ellint_rf_imp(x, y, z) - t * ellint_rd_imp(x, y, z) / 3);
               if(m != 0)
                  result += m * ellint_e_imp(k);
            }
            return invert ? -result : result;
        }

        private static double ellint_e_imp(double k)
        {

            if (Math.Abs(k) > 1)
            {
               throw new Exception(string.Format("Got k = {0:G}, function requires |k| <= 1", k));
            }
            if (Math.Abs(k) == 1)
            {
                return 1.0;
            }

            double x = 0;
            double t = k * k;
            double y = 1 - t;
            double z = 1;
            double value = ellint_rf_imp(x, y, z) - t * ellint_rd_imp(x, y, z) / 3;

            return value;
        }

        private static double ellint_pi_imp(double v, double phi, double k, double vc)
        {
            // Note vc = 1-v presumably without cancellation error.
            double value, x, y, z, p, t;

            if (Math.Abs(k) > 1)
            {
               throw new Exception(string.Format("Got k = {0:G}, function requires |k| <= 1", k));
            }

            double sphi = Math.Sin(Math.Abs(phi));

            if(v > 1 / (sphi * sphi))
            {
                // Complex result is a domain error:
               throw new Exception(string.Format("Got v = {0:G}, but result is complex for v > 1 / Math.Sin^2(phi)", v));
            }

            // Special cases first:
            if(v == 0)
            {
               // A&S 17.7.18 & 19
               return (k == 0) ? phi : ellint_f_imp(phi, k);
            }
            if(phi == Math.PI / 2)
            {
               return ellint_pi_imp(v, k, vc);
            }
            if(k == 0)
            {
               // A&S 17.7.20:
               if(v < 1)
               {
                  double vcr = Math.Sqrt(vc);
                  return Math.Atan(vcr * Math.Tan(phi)) / vcr;
               }
               else if(v == 1)
               {
                  return Math.Tan(phi);
               }
               else
               {
                  // v > 1:
                  double vcr = Math.Sqrt(-vc);
                  double arg = vcr * Math.Tan(phi);
                  return (log1p(arg) - log1p(-arg)) / (2 * vcr);
               }
            }

            if(v < 0)
            {
               //
               // If we don't shift to 0 <= v <= 1 we get
               // cancellation errors later on.  Use
               // A&S 17.7.15/16 to shift to v > 0:
               //
               double k2 = k * k;
               double N = (k2 - v) / (1 - v);
               double Nm1 = (1 - k2) / (1 - v);
               double p2 = Math.Sqrt(-v * (k2 - v) / (1 - v));
               double delta = Math.Sqrt(1 - k2 * sphi * sphi);
               double result = ellint_pi_imp(N, phi, k, Nm1);

               result *= Math.Sqrt(Nm1 * (1 - k2 / N));
               result += ellint_f_imp(phi, k) * k2 / p2;
               result += Math.Atan((p2/2) * Math.Sin(2 * phi) / delta);
               result /= Math.Sqrt((1 - v) * (1 - k2 / v));
               return result;
            }

            if(Math.Abs(phi) > 1 / epsilon)
            {
               if(v > 1) throw new Exception(string.Format("Got v = {0:G}, but this is only supported for 0 <= phi <= pi/2", v));
               value = 2 * Math.Abs(phi) * ellint_pi_imp(v, k, vc) / Math.PI;
            }
            else
            {
               double rphi = fmod(Math.Abs(phi), Math.PI / 2);
               double m = Math.Floor((2 * Math.Abs(phi)) / Math.PI);
               int sign = 1;
               if(fmod(m, 2.0) > 0.5)
               {
                  m += 1;
                  sign = -1;
                  rphi = Math.PI / 2 - rphi;
               }

               double sinp = Math.Sin(rphi);
               double cosp = Math.Cos(rphi);
               x = cosp * cosp;
               t = sinp * sinp;
               y = 1 - k * k * t;
               z = 1;
               if(v * t < 0.5)
                   p = 1 - v * t;
               else
                   p = x + vc * t;
               value = sign * sinp * (ellint_rf_imp(x, y, z) + v * t * ellint_rj_imp(x, y, z, p) / 3);
               if((m > 0) && (vc > 0))
                 value += m * ellint_pi_imp(v, k, vc);
            }

            if (phi < 0)
            {
                value = -value;    // odd function
            }
            return value;
        }

        private static double ellint_pi_imp(double v, double k, double vc)
        {
            if (Math.Abs(k) >= 1)
            {
               throw new Exception(string.Format("Got k = {0:G}, function requires |k| <= 1", k));
            }
            if(vc <= 0)
            {
               // Result is complex:
               throw new Exception(string.Format("Got v = {0:G}, function requires v < 1", v));
            }

            if(v == 0)
            {
               return (k == 0) ? Math.PI / 2 : ellint_k_imp(k);
            }

            if(v < 0)
            {
               double k2 = k * k;
               double N = (k2 - v) / (1 - v);
               double Nm1 = (1 - k2) / (1 - v);
               double p2 = Math.Sqrt(-v * (k2 - v) / (1 - v));

               double result = ellint_pi_imp(N, k, Nm1);

               result *= Math.Sqrt(Nm1 * (1 - k2 / N));
               result += ellint_k_imp(k) * k2 / p2;
               result /= Math.Sqrt((1 - v) * (1 - k2 / v));
               return result;
            }

            double x = 0;
            double y = 1 - k * k;
            double z = 1;
            double p = vc;
            double value = ellint_rf_imp(x, y, z) + v * ellint_rj_imp(x, y, z, p) / 3;

            return value;
        }

        private static void check_series_iterations(ulong k)
        {
            if (k >= max_series_iterations) throw new Exception("Series did not converge");
        }


    }
}
