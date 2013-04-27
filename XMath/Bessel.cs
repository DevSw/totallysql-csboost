using System;
using System.Collections.Generic;
using System.Text;

namespace CSBoost
{
    partial class XMath
    {

        public static double cyl_bessel_i(double v, double x)
        {
            return cyl_bessel_i_imp(v, x);
        }

        public static double cyl_bessel_j(double v, double x)
        {
            return cyl_bessel_j_imp_maybe_int(v, x);
        }

        public static double cyl_bessel_k(double v, double x)
        {
            return cyl_bessel_k_imp_maybe_int(v, x);
        }

        public static double sph_bessel(uint v, double x)
        {
            return sph_bessel_j_imp(v, x);
        }

        public static double cyl_neumann(double v, double x)
        {
            return cyl_neumann_imp_maybe_int(v, x);
        }

        public static double sph_neumann(uint v, double x)
        {
            return sph_neumann_imp(v, x);
        }

        #region implementation

        class bessel_j_small_z_series_term : series<double>
        {
            uint N;
            double v;
            double mult;
            double term;

            public bessel_j_small_z_series_term(double v_, double x)
            {
                N = 0;
                v = v_;
                mult = x / 2;
                term = Math.Pow(mult, v) / gamma(v+1);
                mult *= -mult;
            }

            public override double  next()
            {
                double r = term;
                ++N;
                term *= mult / (N * (N + v));
                return r;
            }
        }

        class sph_bessel_j_small_z_series_term : series<double>
        {
            uint N;
            uint v;
            double mult;
            double term;

            public sph_bessel_j_small_z_series_term(uint v_, double x)
            {
                N = 0;
                v = v_;
                mult = x / 2;
                term = Math.Pow(mult, (double)v) / gamma(v + 1 + 0.5);
                mult *= -mult;
            }

            public override double  next()
            {
                double r = term;
                ++N;
                term *= mult / (N * (N + v + 0.5));
                return r;
            }
        }

        private static double bessel_j_small_z_series(double v, double x)
        {
           bessel_j_small_z_series_term s = new bessel_j_small_z_series_term(v, x);
           int max_iter = XMath.max_series_iterations;
           double result = sum_series(s, XMath.epsilon, ref max_iter, 0);
           return result;
        }

        private static double sph_bessel_j_small_z_series(uint v, double x)
        {
           sph_bessel_j_small_z_series_term s = new sph_bessel_j_small_z_series_term(v, x);
           int max_iter = XMath.max_series_iterations;
           double result = sum_series(s, XMath.epsilon, ref max_iter, 0);
           return result * Math.Sqrt(Math.PI / 4);
        }

        private static double cyl_bessel_j_imp_non_int(double v, double x)
        {
           if(x < 0)
           {
              // better have integer v:
              if(Math.Floor(v) == v)
              {
                 double r = cyl_bessel_j_imp_non_int(v, -x);
                 if(((int)Math.Round(v) & 1) > 0) r = -r;
                 return r;
              }
              else throw new ArgumentException("Cylindrical Bessel_J: x must be >=0 if v is non-integral");
           }
           if(x == 0)
           {
               if (v == 0) return 1;
               if (v > 0) return 0;
               throw new ArgumentException("Cylindrical Bessel_J: v must be >=0 or a negative integer");
           }   
           
           if((v >= 0) && ((x < 1) || (v > x * x / 4) || (x < 5)))
           {
              return bessel_j_small_z_series(v, x);
           }
           
           double j = 0, y = 0;
           bessel_jy(v, x, out j, out y, 1);    //1 = need_j
           return j;
        }

        private static double cyl_bessel_j_imp_maybe_int(double v, double x)
        {
           if((Math.Abs(v) < 200) && (Math.Floor(v) == v))
           {
              return bessel_jn((int)v, x);
           }
           return cyl_bessel_j_imp_non_int(v, x);
        }

        private static double cyl_bessel_j_imp_int(int v, double x)
        {
           return bessel_jn(v, x);
        }

        private static double sph_bessel_j_imp(uint n, double x)
        {
           if(x < 0) throw new ArgumentException("Spherical Bessel_J: x must be >= 0");
           if(n == 0) return sinc_pi(x);
           if(x < 1) return sph_bessel_j_small_z_series(n, x);
           return Math.Sqrt(Math.PI / (2 * x)) * cyl_bessel_j_imp_non_int(n + 0.5, x);
        }

        private static double cyl_bessel_i_imp(double v, double x)
        {
            if(x < 0)
            {
                if(Math.Floor(v) == v)
                {
                    double r = cyl_bessel_i_imp(v, -x);
                    if(((int)Math.Round(v) & 1) > 0) r = -r;
                    return r;
                }
                else throw new ArgumentException("Cylindrical Bessel_I: x must be >=0 if v is non-integral");
            }
            if(x == 0) return (v == 0) ? 1 : 0;
            if(v == 0.5)
            {
                if(x >= log_max_value)
                {
                    double e = Math.Exp(x / 2);
                    return e * (e / Math.Sqrt(2 * x * Math.PI));
                }
                return Math.Sqrt(2 / (x * Math.PI)) * Math.Sinh(x);
            }
            if(v == 0) return bessel_i0(x);
            if(v == 1) return bessel_i1(x);
            double I, K;
            bessel_ik(v, x, out I, out K, 1);       //1 == need_i
            return I;
        }

        private static double cyl_bessel_k_imp_non_int(double v, double x)
        {
           if(x <= 0) throw new ArgumentException("Cylindrical Bessel_K: x must be > 0");
           double I, K;
           bessel_ik(v, x, out I, out K, 2);  //2 == need_k
           return K;
        }

        private static double cyl_bessel_k_imp_maybe_int(double v, double x)
        {
           if((Math.Floor(v) == v)) return bessel_kn((int)v, x);
           return cyl_bessel_k_imp_non_int(v, x);
        }

        private static double cyl_bessel_k_imp_int(int v, double x)
        {
           return bessel_kn(v, x);
        }

        private static double cyl_neumann_imp_non_int(double v, double x)
        {
           if(x <= 0) throw new ArgumentException("Cylindrical Neumann: x must be > 0");
           double j = 0, y = 0;
           bessel_jy(v, x, out j, out y, 2);        //2 == need_y
           return y;
        }

        private static double cyl_neumann_imp_maybe_int(double v, double x)
        {
            double r;
           if(Math.Floor(v) == v)
           {
              if(Math.Abs(x) > 304 && Math.Abs(x) > 5 * Math.Abs(v))     //304 = asymptotic_bessel_y_limit for double precision
              {
                 r = asymptotic_bessel_y_large_x_2(Math.Abs(v), x);
                 if((v < 0) && (((int)v & 1)) > 0) r = -r;
              }
              else r = bessel_yn((int)v, x);
           }
           else r = cyl_neumann_imp_non_int(v, x);
           return r;
        }

        private static double cyl_neumann_imp_int(int v, double x)
        {
           if(Math.Abs(x) > 304 && Math.Abs(x) > 5 * Math.Abs(v))        //304 = asymptotic_bessel_y_limit for double precision
           {
              double r = asymptotic_bessel_y_large_x_2(Math.Abs(v), x);
              if((v < 0) && (v & 1) > 0) r = -r;
              return r;
           }
           else return bessel_yn(v, x);
        }

        private static double sph_neumann_imp(uint v, double x)
        {
            if(x < 0) throw new ArgumentException("Spherical Neumann: x must be >= 0");
            if(x < 2 * double.Epsilon) throw new OverflowException();
            double result = cyl_neumann_imp_non_int(v + 0.5, x);
            double tx = Math.Sqrt(Math.PI / (2 * x));
            if((tx > 1) && (double.MaxValue / tx < result)) throw new OverflowException();
            return result * tx;
        }

        private static bool hankel_PQ(double v, double x, out double p, out double q)
        {
            double tolerance = 2 * XMath.epsilon;
            p = 1;
            q = 0;
            double k = 1;
            double z8 = 8 * x;
            double sq = 1;
            double mu = 4 * v * v;
            double term = 1;
            bool ok = true;
            do
            {
                term *= (mu - sq * sq) / (k * z8);
                q += term;
                k += 1;
                sq += 2;
                double mult = (sq * sq - mu) / (k * z8);
                ok = Math.Abs(mult) < 0.5;
                term *= mult;
                p += term;
                k += 1;
                sq += 2;
            }
            while((Math.Abs(term) > tolerance * p) && ok);
            return ok;
        }

        // Calculate Y(v, x) and Y(v+1, x) by Temme's method, see
        // Temme, Journal of Computational Physics, vol 21, 343 (1976)
        private static int temme_jy(double v, double x, out double Y, out double Y1)
        {
            double g, h, p, q, f, coef, sum, sum1, tolerance;
            double a, d, e, sigma;
            ulong k;

            double gp = tgamma1pm1(v);
            double gm = tgamma1pm1(-v);
            double spv = sin_pi(v);
            double spv2 = sin_pi(v/2);
            double xp = Math.Pow(x/2, v);

            a = Math.Log(x / 2);
            sigma = -a * v;
            d = Math.Abs(sigma) < epsilon ? 1.0 : Math.Sinh(sigma) / sigma;
            e = Math.Abs(v) < epsilon ? v * Math.PI * Math.PI / 2.0 : 2.0 * spv2 * spv2 / v;
            double g1 = (v == 0) ? -euler : (gp - gm) / ((1.0 + gp) * (1.0 + gm) * 2.0 * v);
            double g2 = (2 + gp + gm) / ((1 + gp) * (1 + gm) * 2);
            double vspv = (Math.Abs(v) < epsilon) ? 1/Math.PI : v / spv;
            f = (g1 * Math.Cosh(sigma) - g2 * a * d) * 2 * vspv;

            p = vspv / (xp * (1 + gm));
            q = vspv * xp / (1 + gp);

            g = f + e * q;
            h = p;
            coef = 1;
            sum = coef * g;
            sum1 = coef * h;

            double v2 = v * v;
            double coef_mult = -x * x / 4;

            // series summation
            tolerance = epsilon;
            for (k = 1; k < max_series_iterations; k++)
            {
                f = (k * f + p + q) / (k*k - v2);
                p /= k - v;
                q /= k + v;
                g = f + e * q;
                h = p - k * g;
                coef *= coef_mult / k;
                sum += coef * g;
                sum1 += coef * h;
                if (Math.Abs(coef * g) < Math.Abs(sum) * tolerance) break; 
            }
            if (k == max_series_iterations) throw new Exception();
            Y = -sum;
            Y1 = -2 * sum1 / x;

            return 0;
        }

        // Evaluate continued fraction fv = J_(v+1) / J_v, see
        // Abramowitz and Stegun, Handbook of Mathematical Functions, 1972, 9.1.73
        private static int CF1_jy(double v, double x, out double fv, out int sign)
        {
            double C, D, f, a, b, delta, tiny, tolerance;
            ulong k;
            int s = 1;

            tolerance = 2 * epsilon;
            tiny = Math.Sqrt(double.Epsilon);
            C = f = tiny;                           // b0 = 0, replace with tiny
            D = 0;
            for (k = 1; k < max_series_iterations * 100; k++)
            {
                a = -1;
                b = 2 * (v + k) / x;
                C = b + a / C;
                D = b + a * D;
                if (C == 0) { C = tiny; }
                if (D == 0) { D = tiny; }
                D = 1 / D;
                delta = C * D;
                f *= delta;
                if (D < 0) { s = -s; }
                if (Math.Abs(delta - 1) < tolerance) 
                { break; }
            }
            if (k == max_series_iterations * 100) throw new Exception();
            fv = -f;
            sign = s;                              // sign of denominator

            return 0;
        }

        private static int CF2_jy(double v, double x, out double p, out double q)
        {
            double Cr, Ci, Dr, Di, fr, fi, a, br, bi, delta_r, delta_i, temp;
            double tiny;
            ulong k;

            // modified Lentz's method, complex numbers involved, see
            // Lentz, Applied Optics, vol 15, 668 (1976)
            double tolerance = 2 * epsilon;
            tiny = Math.Sqrt(double.Epsilon);
            Cr = fr = -0.5 / x;
            Ci = fi = 1;
            //Dr = Di = 0;
            double v2 = v * v;
            a = (0.25 - v2) / x; // Note complex this one time only!
            br = 2 * x;
            bi = 2;
            temp = Cr * Cr + 1;
            Ci = bi + a * Cr / temp;
            Cr = br + a / temp;
            Dr = br;
            Di = bi;
            if (Math.Abs(Cr) + Math.Abs(Ci) < tiny) { Cr = tiny; }
            if (Math.Abs(Dr) + Math.Abs(Di) < tiny) { Dr = tiny; }
            temp = Dr * Dr + Di * Di;
            Dr = Dr / temp;
            Di = -Di / temp;
            delta_r = Cr * Dr - Ci * Di;
            delta_i = Ci * Dr + Cr * Di;
            temp = fr;
            fr = temp * delta_r - fi * delta_i;
            fi = temp * delta_i + fi * delta_r;
            for (k = 2; k < max_series_iterations; k++)
            {
                a = k - 0.5;
                a *= a;
                a -= v2;
                bi += 2;
                temp = Cr * Cr + Ci * Ci;
                Cr = br + a * Cr / temp;
                Ci = bi - a * Ci / temp;
                Dr = br + a * Dr;
                Di = bi + a * Di;
                if (Math.Abs(Cr) + Math.Abs(Ci) < tiny) { Cr = tiny; }
                if (Math.Abs(Dr) + Math.Abs(Di) < tiny) { Dr = tiny; }
                temp = Dr * Dr + Di * Di;
                Dr = Dr / temp;
                Di = -Di / temp;
                delta_r = Cr * Dr - Ci * Di;
                delta_i = Ci * Dr + Cr * Di;
                temp = fr;
                fr = temp * delta_r - fi * delta_i;
                fi = temp * delta_i + fi * delta_r;
                if (Math.Abs(delta_r - 1) + Math.Abs(delta_i) < tolerance) break;
            }
            if(k == max_series_iterations) throw new Exception();
            p = fr;
            q = fi;

            return 0;
        }

        const int need_j = 1, need_y = 2;

        // Compute J(v, x) and Y(v, x) simultaneously by Steed's method, see
        // Barnett et al, Computer Physics Communications, vol 8, 377 (1974)
        private static int bessel_jy(double v, double x, out double J, out double Y, int kind)
        {
            double u, Jv, Ju, Yv, Yv1, Yu, Yu1 = 0, fv, fu;
            double W, p, q, gamma, current, prev, next;
            bool reflect = false;
            uint n, k;
            int s;

            if (v < 0)
            {
                reflect = true;
                v = -v;                             // v is non-negative from here
                kind = need_j|need_y;               // need both for reflection formula
            }
            n = (uint)Math.Round(v);
            u = v - n;                              // -1/2 <= u < 1/2

            if (x == 0) throw new OverflowException();

            // x is positive until reflection
            W = 2.0/ (x * Math.PI);               // Wronskian
            if((x > 8) && (x < 1000) && hankel_PQ(v, x, out p, out q))
            {
               //
               // Hankel approximation: note that this method works best when x 
               // is large, but in that case we end up calculating sines and cosines
               // of large values, with horrendous resulting accuracy.  It is fast though
               // when it works....
               //
               double chi = x - fmod(v / 2 + 0.25, 2.0) * Math.PI;
               double sc = sin(chi);
               double cc = cos(chi);
               chi = Math.Sqrt(2 / (Math.PI * x));
               Yv = chi * (p * sc + q * cc);
               Jv = chi * (p * cc - q * sc);
            }
            else if (x <= 2)                           // x in (0, 2]
            {
                if(temme_jy(u, x, out Yu, out Yu1) > 0)             // Temme series
                {
                   // domain error:
                   J = Y = Yu;
                   return 1;
                }
                prev = Yu;
                current = Yu1;
                for (k = 1; k <= n; k++)            // forward recurrence for Y
                {
                    next = 2 * (u + k) * current / x - prev;
                    prev = current;
                    current = next;
                }
                Yv = prev;
                Yv1 = current;
                if((kind & need_j) > 0)
                {
                  CF1_jy(v, x, out fv, out s);                 // continued fraction CF1_jy
                  Jv = W / (Yv * fv - Yv1);           // Wronskian relation
                }
                else
                   Jv = double.NaN; // any value will do, we're not using it.
            }
            else                                    // x in (2, \infty)
            {
                // Get Y(u, x):
                double lim;
                switch(kind)
                {
                case need_j:
                   lim = 33 * Math.Max(3.0, v*v); //asymptotic_bessel_j_limit
                   break;
                case need_y:
                   lim = 304;                     //asymptotic_bessel_y_limit
                   break;
                default:
                   lim = Math.Max(304, 33 * Math.Max(3.0, v*v));
                   break;
                }
                if(x > lim)
                {
                   if((kind & need_y) > 0)
                   {
                      Yu = asymptotic_bessel_y_large_x_2(u, x);
                      Yu1 = asymptotic_bessel_y_large_x_2(u + 1.0, x);
                   }
                   else Yu = double.NaN; // any value will do, we're not using it.
                   if((kind & need_j) > 0) Jv = asymptotic_bessel_j_large_x_2(v, x);
                   else Jv = double.NaN; // any value will do, we're not using it.
                }
                else
                {
                   CF1_jy(v, x, out fv, out s);
                   // tiny initial value to prevent overflow
                   double init = Math.Sqrt(double.Epsilon);
                   prev = fv * s * init;
                   current = s * init;
                   for (k = n; k > 0; k--)             // backward recurrence for J
                   {
                       next = 2 * (u + k) * current / x - prev;
                       prev = current;
                       current = next;
                   }
                   double ratio = (s * init) / current;     // scaling ratio
                   // can also call CF1_jy() to get fu, not much difference in precision
                   fu = prev / current;
                   CF2_jy(u, x, out p, out q);                  // continued fraction CF2_jy
                   double t = u / x - fu;                   // t = J'/J
                   gamma = (p - t) / q;
                   Ju = sign(current) * Math.Sqrt(W / (q + gamma * (p - t)));

                   Jv = Ju * ratio;                    // normalization

                   Yu = gamma * Ju;
                   Yu1 = Yu * (u/x - p - q/gamma);
                }
                if((kind & need_y) > 0)
                {
                   // compute Y:
                   prev = Yu;
                   current = Yu1;
                   for (k = 1; k <= n; k++)            // forward recurrence for Y
                   {
                       next = 2 * (u + k) * current / x - prev;
                       prev = current;
                       current = next;
                   }
                   Yv = prev;
                }
                else
                   Yv = double.NaN; // any value will do, we're not using it.
            }

            if (reflect)
            {
                double z = (u + n % 2);
                J = cos_pi(z) * Jv - sin_pi(z) * Yv;     // reflection formula
                Y = sin_pi(z) * Jv + cos_pi(z) * Yv;
            }
            else
            {
                J = Jv;
                Y = Yv;
            }

            return 0;
        }

        private static double asymptotic_bessel_j_large_x_P(double v, double x)
        {
           // A&S 9.2.9
           double s = 1;
           double mu = 4 * v * v;
           double ez2 = 8 * x;
           ez2 *= ez2;
           s -= (mu-1) * (mu-9) / (2 * ez2);
           s += (mu-1) * (mu-9) * (mu-25) * (mu - 49) / (24 * ez2 * ez2);
           return s;
        }

        private static double asymptotic_bessel_j_large_x_Q(double v, double x)
        {
           // A&S 9.2.10
           double s = 0;
           double mu = 4 * v * v;
           double ez = 8*x;
           s += (mu-1) / ez;
           s -= (mu-1) * (mu-9) * (mu-25) / (6 * ez*ez*ez);
           return s;
        }

        private static double asymptotic_bessel_j_large_x(double v, double x)
        {
           // 
           // See http://functions.wolfram.com/BesselAiryStruveFunctions/BesselJ/06/02/02/0001/
           //
           // Also A&S 9.2.5
           //
           double chi = Math.Abs(x) - Math.PI * (2 * v + 1) / 4;
           return Math.Sqrt(2 / (Math.PI * x))
              * (asymptotic_bessel_j_large_x_P(v, x) * cos(chi) 
                 - asymptotic_bessel_j_large_x_Q(v, x) * sin(chi));
        }

        private static double asymptotic_bessel_y_large_x(double v, double x)
        {
           // 
           // See http://functions.wolfram.com/BesselAiryStruveFunctions/BesselJ/06/02/02/0001/
           //
           // Also A&S 9.2.5
           //
           double chi = Math.Abs(x) - Math.PI * (2 * v + 1) / 4;
           return Math.Sqrt(2 / (Math.PI * x))
              * (asymptotic_bessel_j_large_x_P(v, x) * sin(chi) 
                 - asymptotic_bessel_j_large_x_Q(v, x) * cos(chi));
        }

        private static double asymptotic_bessel_amplitude(double v, double x)
        {
           // Calculate the amplitude of J(v, x) and Y(v, x) for large
           // x: see A&S 9.2.28.
           double s = 1;
           double mu = 4 * v * v;
           double txq = 2 * x;
           txq *= txq;

           s += (mu - 1) / (2 * txq);
           s += 3 * (mu - 1) * (mu - 9) / (txq * txq * 8);
           s += 15 * (mu - 1) * (mu - 9) * (mu - 25) / (txq * txq * txq * 8 * 6);

           return Math.Sqrt(s * 2 / (Math.PI * x));
        }

        private static double asymptotic_bessel_phase_mx(double v, double x)
        {
           //
           // Calculate the phase of J(v, x) and Y(v, x) for large x.
           // See A&S 9.2.29.
           // Note that the result returned is the phase less x.
           //
           double mu = 4 * v * v;
           double denom = 4 * x;
           double denom_mult = denom * denom;

           double s = -Math.PI * (v / 2 + 0.25f);
           s += (mu - 1) / (2 * denom);
           denom *= denom_mult;
           s += (mu - 1) * (mu - 25) / (6 * denom);
           denom *= denom_mult;
           s += (mu - 1) * (mu * mu - 114 * mu + 1073) / (5 * denom);
           denom *= denom_mult;
           s += (mu - 1) * (5 * mu * mu * mu - 1535 * mu * mu + 54703 * mu - 375733) / (14 * denom);
           return s;
        }

        private static double asymptotic_bessel_y_large_x_2(double v, double x)
        {
            // See A&S 9.2.19.
            // Get the phase and amplitude:
            double ampl = asymptotic_bessel_amplitude(v, x);
            double phase = asymptotic_bessel_phase_mx(v, x);
            //
            // Calculate the sine of the phase, using:
            // sin(x+p) = sin(x)cos(p) + cos(x)sin(p)
            //
            double sin_phase = sin(phase) * cos(x) + cos(phase) * sin(x);
            return sin_phase * ampl;
        }

        private static double asymptotic_bessel_j_large_x_2(double v, double x)
        {
           // See A&S 9.2.19.
           // Get the phase and amplitude:
           double ampl = asymptotic_bessel_amplitude(v, x);
           double phase = asymptotic_bessel_phase_mx(v, x);
           //
           // Calculate the sine of the phase, using:
           // cos(x+p) = cos(x)cos(p) - sin(x)sin(p)
           //
           double sin_phase = cos(phase) * cos(x) - sin(phase) * sin(x);
           return sin_phase * ampl;
        }
        
        // Calculate K(v, x) and K(v+1, x) by method analogous to
        // Temme, Journal of Computational Physics, vol 21, 343 (1976)
        private static int temme_ik(double v, double x, out double K, out double K1)
        {
            double f, h, p, q, coef, sum, sum1, tolerance;
            double a, b, c, d, sigma, gamma1, gamma2;
            ulong k;

            double gp = tgamma1pm1(v);
            double gm = tgamma1pm1(-v);

            a = Math.Log(x / 2);
            b = Math.Exp(v * a);
            sigma = -a * v;
            c = Math.Abs(v) < epsilon ? 1.0 : sin_pi(v) / (v * Math.PI);
            d = Math.Abs(sigma) < epsilon ? 1.0 : Math.Sinh(sigma) / sigma;
            gamma1 = Math.Abs(v) < epsilon ? -euler : (0.5 / v) * (gp - gm) * c;
            gamma2 = (2 + gp + gm) * c / 2;

            // initial values
            p = (gp + 1) / (2 * b);
            q = (1 + gm) * b / 2;
            f = (Math.Cosh(sigma) * gamma1 + d * (-a) * gamma2) / c;
            h = p;
            coef = 1;
            sum = coef * f;
            sum1 = coef * h;

            // series summation
            tolerance = epsilon;
            for (k = 1; k < max_series_iterations; k++)
            {
                f = (k * f + p + q) / (k*k - v*v);
                p /= k - v;
                q /= k + v;
                h = p - k * f;
                coef *= x * x / (4 * k);
                sum += coef * f;
                sum1 += coef * h;
                if (Math.Abs(coef * f) < Math.Abs(sum) * tolerance) 
                { 
                   break; 
                }
            }
            if(k == max_series_iterations) throw new Exception();

            K = sum;
            K1 = 2 * sum1 / x;

            return 0;
        }

        // Evaluate continued fraction fv = I_(v+1) / I_v, derived from
        // Abramowitz and Stegun, Handbook of Mathematical Functions, 1972, 9.1.73
        private static int CF1_ik(double v, double x, out double fv)
        {
            double C, D, f, a, b, delta, tiny, tolerance;
            ulong k;

            // modified Lentz's method, see
            // Lentz, Applied Optics, vol 15, 668 (1976)
            tolerance = 2 * epsilon;
            tiny = Math.Sqrt(double.Epsilon);
            C = f = tiny;                           // b0 = 0, replace with tiny
            D = 0;
            for (k = 1; k < max_series_iterations; k++)
            {
                a = 1;
                b = 2 * (v + k) / x;
                C = b + a / C;
                D = b + a * D;
                if (C == 0) { C = tiny; }
                if (D == 0) { D = tiny; }
                D = 1 / D;
                delta = C * D;
                f *= delta;
                if (Math.Abs(delta - 1) <= tolerance) 
                { 
                   break; 
                }
            }
            if(k == max_series_iterations) throw new Exception();
            fv = f;
            return 0;
        }

        // Calculate K(v, x) and K(v+1, x) by evaluating continued fraction
        // z1 / z0 = U(v+1.5, 2v+1, 2x) / U(v+0.5, 2v+1, 2x), see
        // Thompson and Barnett, Computer Physics Communications, vol 47, 245 (1987)
        private static int CF2_ik(double v, double x, out double Kv, out double Kv1)
        {
            double S, C, Q, D, f, a, b, q, delta, tolerance, current, prev;
            ulong k;

            // Steed's algorithm, see Thompson and Barnett,
            // Journal of Computational Physics, vol 64, 490 (1986)
            tolerance = epsilon;
            a = v * v - 0.25f;
            b = 2 * (x + 1);                              // b1
            D = 1 / b;                                    // D1 = 1 / b1
            f = delta = D;                                // f1 = delta1 = D1, coincidence
            prev = 0;                                     // q0
            current = 1;                                  // q1
            Q = C = -a;                                   // Q1 = C1 because q1 = 1
            S = 1 + Q * delta;                            // S1
            for (k = 2; k < max_series_iterations; k++)     // starting from 2
            {
                // continued fraction f = z1 / z0
                a -= 2 * (k - 1);
                b += 2;
                D = 1 / (b + a * D);
                delta *= b * D - 1;
                f += delta;

                // series summation S = 1 + \sum_{n=1}^{\infty} C_n * z_n / z_0
                q = (prev - (b - 2) * current) / a;
                prev = current;
                current = q;                        // forward recurrence for q
                C *= -a / k;
                Q += C * q;
                S += Q * delta;

                // S converges slower than f
                if (Math.Abs(Q * delta) < Math.Abs(S) * tolerance) break; 
            }
            if (k == max_series_iterations) throw new Exception();

            Kv = Math.Sqrt(Math.PI / (2 * x)) * Math.Exp(-x) / S;
            Kv1 = Kv * (0.5 + v + x + (v * v - 0.25) * f) / x;
            return 0;
        }

        const int need_i = 1, need_k = 2;

        // Compute I(v, x) and K(v, x) simultaneously by Temme's method, see
        // Temme, Journal of Computational Physics, vol 19, 324 (1975)
        private static int bessel_ik(double v, double x, out double I, out double K, int kind)
        {
            double u, Iv, Kv, Kv1, Ku, Ku1, fv;
            double W, current, prev, next;
            bool reflect = false;
            uint n, k;

            if (v < 0)
            {
                reflect = true;
                v = -v;                             // v is non-negative from here
                kind |= need_k;
            }
            n = (uint)Math.Round(v);
            u = v - n;                              // -1/2 <= u < 1/2

            if (x < 0) throw new ArgumentException("x must be >= 0");
            if (x == 0)
            {
               Iv = (v == 0) ? 1.0: 0.0;
               if((kind & need_k) > 0) throw new OverflowException();
               Kv = double.NaN; // any value will do

               if(reflect && (kind & need_i) > 0)
               {
                   double z = (u + n % 2);
                   Iv = sin_pi(z);
                   if (Iv == 0) throw new OverflowException();
               }

               I = Iv;
               K = Kv;
               return 0;
            }

            // x is positive until reflection
            W = 1 / x;                                 // Wronskian
            if (x <= 2)                                // x in (0, 2]
            {
                temme_ik(u, x, out Ku, out Ku1);             // Temme series
            }
            else                                       // x in (2, \infty)
            {
                CF2_ik(u, x, out Ku, out Ku1);               // continued fraction CF2_ik
            }
            prev = Ku;
            current = Ku1;
            for (k = 1; k <= n; k++)                   // forward recurrence for K
            {
                next = 2 * (u + k) * current / x + prev;
                prev = current;
                current = next;
            }
            Kv = prev;
            Kv1 = current;
            if((kind & need_i) > 0)
            {
               double lim = (4 * v * v + 10) / (8 * x);
               lim *= lim;
               lim *= lim;
               lim /= 24;
               if((lim < epsilon * 10) && (x > 100))
               {
                  // x is huge compared to v, CF1 may be very slow
                  // to converge so use asymptotic expansion for large
                  // x case instead.  Note that the asymptotic expansion
                  // isn't very accurate - so it's deliberately very hard 
                  // to get here - probably we're going to overflow:
                  Iv = asymptotic_bessel_i_large_x(v, x);
               }
               else
               {
                  CF1_ik(v, x, out fv);                         // continued fraction CF1_ik
                  Iv = W / (Kv * fv + Kv1);                  // Wronskian relation
               }
            }
            else Iv = double.NaN; // any value will do

            if (reflect)
            {
                double z = (u + n % 2);
                I = Iv + (2 / Math.PI) * sin_pi(z) * Kv;   // reflection formula
                K = Kv;
            }
            else
            {
                I = Iv;
                K = Kv;
            }
            return 0;
        }

        private static double asymptotic_bessel_i_large_x(double v, double x)
        {
            double s = 1;
            double mu = 4 * v * v;
            double ex = 8 * x;
            double num = mu - 1;
            double denom = ex;

            s -= num / denom;

            num *= mu - 9;
            denom *= ex * 2;
            s += num / denom;

            num *= mu - 25;
            denom *= ex * 3;
            s -= num / denom;

            // Try and avoid overflow to the last minute:
            double e = Math.Exp(x / 2);

            s = e * (e * s / Math.Sqrt(2 * x * Math.PI));

            return s;
        }

        private static double bessel_jn(int n, double x)
        {
            double value = 0, factor, current, prev, next;

            if (n < 0)
            {
                factor = ((n & 0x1) > 0) ? -1 : 1;  // J_{-n}(z) = (-1)^n J_n(z)
                n = -n;
            }
            else factor = 1;
            //
            // Special cases:
            //
            if (n == 0) return factor * bessel_j0(x);
            if (n == 1) return factor * bessel_j1(x);
            if (x == 0) return 0.0;

            if (Math.Abs(x) > 33 * Math.Max(3.0, n * n)) //asymptotic_bessel_j_limit
                return factor * asymptotic_bessel_j_large_x_2(n, x);

            if (n < Math.Abs(x))                         // forward recurrence
            {
                prev = bessel_j0(x);
                current = bessel_j1(x);
                for (int k = 1; k < n; k++)
                {
                    value = 2 * k * current / x - prev;
                    prev = current;
                    current = value;
                }
            }
            else                                    // backward recurrence
            {
                double fn; int s;                        // fn = J_(n+1) / J_n
                // |x| <= n, fast convergence for continued fraction CF1
                CF1_jy(n, x, out fn, out s);
                // tiny initial value to prevent overflow
                double init = Math.Sqrt(double.Epsilon);
                prev = fn * init;
                current = init;
                for (int k = n; k > 0; k--)
                {
                    next = 2 * k * current / x - prev;
                    prev = current;
                    current = next;
                }
                double ratio = init / current;           // scaling ratio
                value = bessel_j0(x) * ratio;       // normalization
            }
            value *= factor;

            return value;
        }

        private static double bessel_j0(double x)
        {
            double[] P1 = 
            {
                 -4.1298668500990866786e+11,
                 2.7282507878605942706e+10,
                 -6.2140700423540120665e+08,
                 6.6302997904833794242e+06,
                 -3.6629814655107086448e+04,
                 1.0344222815443188943e+02,
                 -1.2117036164593528341e-01
            };
            double[] Q1 = 
            {
                 2.3883787996332290397e+12,
                 2.6328198300859648632e+10,
                 1.3985097372263433271e+08,
                 4.5612696224219938200e+05,
                 9.3614022392337710626e+02,
                 1.0,
                 0.0
            };
            double[] P2 = 
            {
                 -1.8319397969392084011e+03,
                 -1.2254078161378989535e+04,
                 -7.2879702464464618998e+03,
                 1.0341910641583726701e+04,
                 1.1725046279757103576e+04,
                 4.4176707025325087628e+03,
                 7.4321196680624245801e+02,
                 4.8591703355916499363e+01
            };
            double[] Q2 = 
            {
                 -3.5783478026152301072e+05,
                 2.4599102262586308984e+05,
                 -8.4055062591169562211e+04,
                 1.8680990008359188352e+04,
                 -2.9458766545509337327e+03,
                 3.3307310774649071172e+02,
                 -2.5258076240801555057e+01,
                 1.0
            };
            double[] PC = 
            {
                 2.2779090197304684302e+04,
                 4.1345386639580765797e+04,
                 2.1170523380864944322e+04,
                 3.4806486443249270347e+03,
                 1.5376201909008354296e+02,
                 8.8961548424210455236e-01
            };
            double[] QC = 
            {
                 2.2779090197304684318e+04,
                 4.1370412495510416640e+04,
                 2.1215350561880115730e+04,
                 3.5028735138235608207e+03,
                 1.5711159858080893649e+02,
                 1.0
            };
            double[] PS = 
            {
                -8.9226600200800094098e+01,
                -1.8591953644342993800e+02,
                -1.1183429920482737611e+02,
                -2.2300261666214198472e+01,
                -1.2441026745835638459e+00,
                -8.8033303048680751817e-03
            };
            double[] QS = 
            {
                 5.7105024128512061905e+03,
                 1.1951131543434613647e+04,
                 7.2642780169211018836e+03,
                 1.4887231232283756582e+03,
                 9.0593769594993125859e+01,
                 1.0
            };

            double x1 = 2.4048255576957727686e+00,
                    x2 = 5.5200781102863106496e+00,
                    x11 = 6.160e+02,
                    x12 = -1.42444230422723137837e-03,
                    x21 = 1.4130e+03,
                    x22 = 5.46860286310649596604e-04;

            double value, factor, r, rc, rs;

            if (x < 0) x = -x;                         // even function
            if (x == 0) return 1.0;
            if (x <= 4)
            {
                double y = x * x;
                r = evaluate_rational(P1, Q1, y);
                factor = (x + x1) * ((x - x11 / 256) - x12);
                value = factor * r;
            }
            else if (x <= 8.0)                  // x in (4, 8]
            {
                double y = 1 - (x * x) / 64;
                r = evaluate_rational(P2, Q2, y);
                factor = (x + x2) * ((x - x21 / 256) - x22);
                value = factor * r;
            }
            else                                // x in (8, \infty)
            {
                double y = 8 / x;
                double y2 = y * y;
                double z = x - 0.25 * Math.PI;
                rc = evaluate_rational(PC, QC, y2);
                rs = evaluate_rational(PS, QS, y2);
                factor = Math.Sqrt(2 / (x * Math.PI));
                value = factor * (rc * cos(z) - y * rs * sin(z));
            }

            return value;
        }

        private static double bessel_j1(double x)
        {
            double[] P1 = 
            {
                -1.4258509801366645672e+11,
                6.6781041261492395835e+09,
                -1.1548696764841276794e+08,
                9.8062904098958257677e+05,
                -4.4615792982775076130e+03,
                1.0650724020080236441e+01,
                -1.0767857011487300348e-02
            };
            double[] Q1 = 
            {
                4.1868604460820175290e+12,
                4.2091902282580133541e+10,
                2.0228375140097033958e+08,
                5.9117614494174794095e+05,
                1.0742272239517380498e+03,
                1.0,
                0.0
            };
            double[] P2 = 
            {
         -1.7527881995806511112e+16,
         1.6608531731299018674e+15,
         -3.6658018905416665164e+13,
         3.5580665670910619166e+11,
         -1.8113931269860667829e+09,
         5.0793266148011179143e+06,
         -7.5023342220781607561e+03,
         4.6179191852758252278e+00
    };
            double[] Q2 = {
         1.7253905888447681194e+18,
         1.7128800897135812012e+16,
         8.4899346165481429307e+13,
         2.7622777286244082666e+11,
         6.4872502899596389593e+08,
         1.1267125065029138050e+06,
         1.3886978985861357615e+03,
         1.0
    };
            double[] PC = {
        -4.4357578167941278571e+06,
        -9.9422465050776411957e+06,
        -6.6033732483649391093e+06,
        -1.5235293511811373833e+06,
        -1.0982405543459346727e+05,
        -1.6116166443246101165e+03,
        0.0
    };
            double[] QC = {
        -4.4357578167941278568e+06,
        -9.9341243899345856590e+06,
        -6.5853394797230870728e+06,
        -1.5118095066341608816e+06,
        -1.0726385991103820119e+05,
        -1.4550094401904961825e+03,
        1.0
    };
            double[] PS = {
         3.3220913409857223519e+04,
         8.5145160675335701966e+04,
         6.6178836581270835179e+04,
         1.8494262873223866797e+04,
         1.7063754290207680021e+03,
         3.5265133846636032186e+01,
         0.0
    };
            double[] QS = {
         7.0871281941028743574e+05,
         1.8194580422439972989e+06,
         1.4194606696037208929e+06,
         4.0029443582266975117e+05,
         3.7890229745772202641e+04,
         8.6383677696049909675e+02,
         1.0
    };
            double x1 = 3.8317059702075123156e+00,
                    x2 = 7.0155866698156187535e+00,
                    x11 = 9.810e+02,
                    x12 = -3.2527979248768438556e-04,
                    x21 = 1.7960e+03,
                    x22 = -3.8330184381246462950e-05;

            double value, factor, r, rc, rs, w;

            w = Math.Abs(x);
            if (x == 0) return 0.0;
            if (w <= 4)                       // w in (0, 4]
            {
                double y = x * x;
                r = evaluate_rational(P1, Q1, y);
                factor = w * (w + x1) * ((w - x11 / 256) - x12);
                value = factor * r;
            }
            else if (w <= 8)                  // w in (4, 8]
            {
                double y = x * x;
                r = evaluate_rational(P2, Q2, y);
                factor = w * (w + x2) * ((w - x21 / 256) - x22);
                value = factor * r;
            }
            else                                // w in (8, \infty)
            {
                double y = 8 / w;
                double y2 = y * y;
                double z = w - 0.75 * Math.PI;
                rc = evaluate_rational(PC, QC, y2);
                rs = evaluate_rational(PS, QS, y2);
                factor = Math.Sqrt(2 / (w * Math.PI));
                value = factor * (rc * cos(z) - y * rs * sin(z));
            }

            if (x < 0)
            {
                value *= -1;                 // odd function
            }
            return value;
        }

        private static double bessel_i0(double x)
        {
            double[] P1 = 
            {
                -2.2335582639474375249e+15,
                -5.5050369673018427753e+14,
                -3.2940087627407749166e+13,
                -8.4925101247114157499e+11,
                -1.1912746104985237192e+10,
                -1.0313066708737980747e+08,
                -5.9545626019847898221e+05,
                -2.4125195876041896775e+03,
                -7.0935347449210549190e+00,
                -1.5453977791786851041e-02,
                -2.5172644670688975051e-05,
                -3.0517226450451067446e-08,
                -2.6843448573468483278e-11,
                -1.5982226675653184646e-14,
                -5.2487866627945699800e-18,
            };
            double[] Q1 = 
            {
                -2.2335582639474375245e+15,
                7.8858692566751002988e+12,
                -1.2207067397808979846e+10,
                1.0377081058062166144e+07,
                -4.8527560179962773045e+03,
                1.0,
            };
            double[] P2 = 
            {
                -2.2210262233306573296e-04,
                1.3067392038106924055e-02,
                -4.4700805721174453923e-01,
                5.5674518371240761397e+00,
                -2.3517945679239481621e+01,
                3.1611322818701131207e+01,
                -9.6090021968656180000e+00,
            };
            double[] Q2 = 
            {
                -5.5194330231005480228e-04,
                3.2547697594819615062e-02,
                -1.1151759188741312645e+00,
                1.3982595353892851542e+01,
                -6.0228002066743340583e+01,
                8.5539563258012929600e+01,
                -3.1446690275135491500e+01,
                1.0,
            };
            double value, factor, r;

            if (x < 0) x = -x;                         // even function
            if (x == 0) return 1;
            if (x <= 15)                        // x in (0, 15]
            {
                double y = x * x;
                value = evaluate_polynomial(P1, y) / evaluate_polynomial(Q1, y);
            }
            else                                // x in (15, \infty)
            {
                double y = 1 / x - 1.0 / 15;
                r = evaluate_polynomial(P2, y) / evaluate_polynomial(Q2, y);
                factor = Math.Exp(x) / Math.Sqrt(x);
                value = factor * r;
            }
            return value;
        }

        private static double bessel_i1(double x)
        {
            double[] P1= 
            {
                -1.4577180278143463643e+15,
                -1.7732037840791591320e+14,
                -6.9876779648010090070e+12,
                -1.3357437682275493024e+11,
                -1.4828267606612366099e+09,
                -1.0588550724769347106e+07,
                -5.1894091982308017540e+04,
                -1.8225946631657315931e+02,
                -4.7207090827310162436e-01,
                -9.1746443287817501309e-04,
                -1.3466829827635152875e-06,
                -1.4831904935994647675e-09,
                -1.1928788903603238754e-12,
                -6.5245515583151902910e-16,
                -1.9705291802535139930e-19,
            };
            double[] Q1= 
            {
                -2.9154360556286927285e+15,
                9.7887501377547640438e+12,
                -1.4386907088588283434e+10,
                1.1594225856856884006e+07,
                -5.1326864679904189920e+03,
                1.0,
            };
            double[] P2= 
            {
                1.4582087408985668208e-05,
                -8.9359825138577646443e-04,
                2.9204895411257790122e-02,
                -3.4198728018058047439e-01,
                1.3960118277609544334e+00,
                -1.9746376087200685843e+00,
                8.5591872901933459000e-01,
                -6.0437159056137599999e-02,
            };
            double[] Q2= 
            {
                3.7510433111922824643e-05,
                -2.2835624489492512649e-03,
                7.4212010813186530069e-02,
                -8.5017476463217924408e-01,
                3.2593714889036996297e+00,
                -3.8806586721556593450e+00,
                1.0,
            };
            double value, factor, r, w;

            w = Math.Abs(x);
            if (x == 0)
            {
                return 0;
            }
            if (w <= 15)                        // w in (0, 15]
            {
                double y = x * x;
                r = evaluate_polynomial(P1, y) / evaluate_polynomial(Q1, y);
                factor = w;
                value = factor * r;
            }
            else                                // w in (15, \infty)
            {
                double y = 1 / w - 1.0 / 15;
                r = evaluate_polynomial(P2, y) / evaluate_polynomial(Q2, y);
                factor = Math.Exp(w) / Math.Sqrt(w);
                value = factor * r;
            }

            if (x < 0)
            {
                value *= -value;                 // odd function
            }
            return value;
        }

        private static double bessel_kn(int n, double x)
        {
            double value, current, prev;

            if (x <= 0) throw new Exception("x must be > 0");

            if (n < 0) n = -n;                             // K_{-n}(z) = K_n(z)
            if (n == 0) value = bessel_k0(x);
            else if (n == 1) value = bessel_k1(x);
            else
            {
                prev = bessel_k0(x);
                current = bessel_k1(x);
                int k = 1;
                do
                {
                    value = 2 * k * current / x + prev;
                    prev = current;
                    current = value;
                    ++k;
                }
                while (k < n);
            }
            return value;
        }

        private static double bessel_k0(double x)
        {
            double[] P1 = 
            {
                 2.4708152720399552679e+03,
                 5.9169059852270512312e+03,
                 4.6850901201934832188e+02,
                 1.1999463724910714109e+01,
                 1.3166052564989571850e-01,
                 5.8599221412826100000e-04
            };
            double[] Q1 = 
            {
                 2.1312714303849120380e+04,
                -2.4994418972832303646e+02,
                 1.0
            };
            double[] P2 = 
            {
                -1.6128136304458193998e+06,
                -3.7333769444840079748e+05,
                -1.7984434409411765813e+04,
                -2.9501657892958843865e+02,
                -1.6414452837299064100e+00
            };
            double[] Q2 = 
            {
                -1.6128136304458193998e+06,
                2.9865713163054025489e+04,
                -2.5064972445877992730e+02,
                1.0
            };
            double[] P3 = 
            {
                 1.1600249425076035558e+02,
                 2.3444738764199315021e+03,
                 1.8321525870183537725e+04,
                 7.1557062783764037541e+04,
                 1.5097646353289914539e+05,
                 1.7398867902565686251e+05,
                 1.0577068948034021957e+05,
                 3.1075408980684392399e+04,
                 3.6832589957340267940e+03,
                 1.1394980557384778174e+02
            };
            double[] Q3 = 
            {
                 9.2556599177304839811e+01,
                 1.8821890840982713696e+03,
                 1.4847228371802360957e+04,
                 5.8824616785857027752e+04,
                 1.2689839587977598727e+05,
                 1.5144644673520157801e+05,
                 9.7418829762268075784e+04,
                 3.1474655750295278825e+04,
                 4.4329628889746408858e+03,
                 2.0013443064949242491e+02,
                 1.0
            };
            double value, factor, r, r1, r2;

            if (x <= 0) throw new Exception("x must be > 0");
            if (x <= 1)                         // x in (0, 1]
            {
                double y = x * x;
                r1 = evaluate_polynomial(P1, y) / evaluate_polynomial(Q1, y);
                r2 = evaluate_polynomial(P2, y) / evaluate_polynomial(Q2, y);
                factor = Math.Log(x);
                value = r1 - factor * r2;
            }
            else                                // x in (1, \infty)
            {
                double y = 1 / x;
                r = evaluate_polynomial(P3, y) / evaluate_polynomial(Q3, y);
                factor = Math.Exp(-x) / Math.Sqrt(x);
                value = factor * r;
            }

            return value;
        }

        private static double bessel_k1(double x)
        {
            double[] P1 = 
            {
                -2.2149374878243304548e+06,
                 7.1938920065420586101e+05,
                 1.7733324035147015630e+05,
                 7.1885382604084798576e+03,
                 9.9991373567429309922e+01,
                 4.8127070456878442310e-01
            };
            double[] Q1 = 
            {
                -2.2149374878243304548e+06,
                 3.7264298672067697862e+04,
                -2.8143915754538725829e+02,
                 1.0
            };
            double[] P2 = 
            {
                 0.0,
                -1.3531161492785421328e+06,
                -1.4758069205414222471e+05,
                -4.5051623763436087023e+03,
                -5.3103913335180275253e+01,
                -2.2795590826955002390e-01
            };
            double[] Q2 = 
            {
                -2.7062322985570842656e+06,
                4.3117653211351080007e+04,
                -3.0507151578787595807e+02,
                1.0
            };
            double[] P3 = 
            {
                 2.2196792496874548962e+00,
                 4.4137176114230414036e+01,
                 3.4122953486801312910e+02,
                 1.3319486433183221990e+03,
                 2.8590657697910288226e+03,
                 3.4540675585544584407e+03,
                 2.3123742209168871550e+03,
                 8.1094256146537402173e+02,
                 1.3182609918569941308e+02,
                 7.5584584631176030810e+00,
                 6.4257745859173138767e-02
            };
            double[] Q3 = 
            {
                 1.7710478032601086579e+00,
                 3.4552228452758912848e+01,
                 2.5951223655579051357e+02,
                 9.6929165726802648634e+02,
                 1.9448440788918006154e+03,
                 2.1181000487171943810e+03,
                 1.2082692316002348638e+03,
                 3.3031020088765390854e+02,
                 3.6001069306861518855e+01,
                 1.0
            };
            double value, factor, r, r1, r2;

            if (x < 0) throw new Exception("x must be < 0");
            if (x <= 1)                         // x in (0, 1]
            {
                double y = x * x;
                r1 = evaluate_polynomial(P1, y) / evaluate_polynomial(Q1, y);
                r2 = evaluate_polynomial(P2, y) / evaluate_polynomial(Q2, y);
                factor = Math.Log(x);
                value = (r1 + factor * r2) / x;
            }
            else                                // x in (1, \infty)
            {
                double y = 1 / x;
                r = evaluate_polynomial(P3, y) / evaluate_polynomial(Q3, y);
                factor = Math.Exp(-x) / Math.Sqrt(x);
                value = factor * r;
            }

            return value;
        }

        private static double bessel_yn(int n, double x)
        {
            double value, factor, current, prev;

            if (x <= 0) throw new ArgumentException("x must be > 0)");

            //
            // Reflection comes first:
            //
            if (n < 0)
            {
                factor = ((n & 0x1) > 0) ? -1 : 1;  // Y_{-n}(z) = (-1)^n Y_n(z)
                n = -n;
            }
            else factor = 1;

            if (n == 0) value = bessel_y0(x);
            else if (n == 1) value = factor * bessel_y1(x);
            else
            {
                prev = bessel_y0(x);
                current = bessel_y1(x);
                int k = 1;
                do
                {
                    value = 2 * k * current / x - prev;
                    prev = current;
                    current = value;
                    ++k;
                }
                while (k < n);
                value *= factor;
            }
            return value;
        }

        private static double bessel_y0(double x)
        {
            double[] P1 = 
            {
                 1.0723538782003176831e+11,
                -8.3716255451260504098e+09,
                 2.0422274357376619816e+08,
                -2.1287548474401797963e+06,
                 1.0102532948020907590e+04,
                -1.8402381979244993524e+01,
            };
            double[] Q1 = 
            {
                 5.8873865738997033405e+11,
                 8.1617187777290363573e+09,
                 5.5662956624278251596e+07,
                 2.3889393209447253406e+05,
                 6.6475986689240190091e+02,
                 1.0,
            };
            double[] P2 = 
            {
                -2.2213976967566192242e+13,
                -5.5107435206722644429e+11,
                 4.3600098638603061642e+10,
                -6.9590439394619619534e+08,
                 4.6905288611678631510e+06,
                -1.4566865832663635920e+04,
                 1.7427031242901594547e+01,
            };
            double[] Q2 = 
            {
                 4.3386146580707264428e+14,
                 5.4266824419412347550e+12,
                 3.4015103849971240096e+10,
                 1.3960202770986831075e+08,
                 4.0669982352539552018e+05,
                 8.3030857612070288823e+02,
                 1.0,
            };
            double[] P3 = 
            {
                -8.0728726905150210443e+15,
                 6.7016641869173237784e+14,
                -1.2829912364088687306e+11,
                -1.9363051266772083678e+11,
                 2.1958827170518100757e+09,
                -1.0085539923498211426e+07,
                 2.1363534169313901632e+04,
                -1.7439661319197499338e+01,
            };
            double[] Q3 = 
            {
                 3.4563724628846457519e+17,
                 3.9272425569640309819e+15,
                 2.2598377924042897629e+13,
                 8.6926121104209825246e+10,
                 2.4727219475672302327e+08,
                 5.3924739209768057030e+05,
                 8.7903362168128450017e+02,
                 1.0,
            };
            double[] PC = 
            {
                 2.2779090197304684302e+04,
                 4.1345386639580765797e+04,
                 2.1170523380864944322e+04,
                 3.4806486443249270347e+03,
                 1.5376201909008354296e+02,
                 8.8961548424210455236e-01,
            };
            double[] QC = 
            {
                 2.2779090197304684318e+04,
                 4.1370412495510416640e+04,
                 2.1215350561880115730e+04,
                 3.5028735138235608207e+03,
                 1.5711159858080893649e+02,
                 1.0,
            };
            double[] PS = 
            {
                -8.9226600200800094098e+01,
                -1.8591953644342993800e+02,
                -1.1183429920482737611e+02,
                -2.2300261666214198472e+01,
                -1.2441026745835638459e+00,
                -8.8033303048680751817e-03,
            };
            double[] QS = 
            {
                 5.7105024128512061905e+03,
                 1.1951131543434613647e+04,
                 7.2642780169211018836e+03,
                 1.4887231232283756582e+03,
                 9.0593769594993125859e+01,
                 1.0,
            };

            double x1 = 8.9357696627916752158e-01,
                    x2 = 3.9576784193148578684e+00,
                    x3 = 7.0860510603017726976e+00,
                    x11 = 2.280e+02,
                    x12 = 2.9519662791675215849e-03,
                    x21 = 1.0130e+03,
                    x22 = 6.4716931485786837568e-04,
                    x31 = 1.8140e+03,
                    x32 = 1.1356030177269762362e-04;

            double value, factor, r, rc, rs;

            if (x <= 0) throw new ArgumentException("x must be > 0");
            if (x <= 3)                       // x in (0, 3]
            {
                double y = x * x;
                double z = 2 * Math.Log(x / x1) * bessel_j0(x) / Math.PI;
                r = evaluate_rational(P1, Q1, y);
                factor = (x + x1) * ((x - x11 / 256) - x12);
                value = z + factor * r;
            }
            else if (x <= 5.5)                  // x in (3, 5.5]
            {
                double y = x * x;
                double z = 2 * Math.Log(x / x2) * bessel_j0(x) / Math.PI;
                r = evaluate_rational(P2, Q2, y);
                factor = (x + x2) * ((x - x21 / 256) - x22);
                value = z + factor * r;
            }
            else if (x <= 8)                  // x in (5.5, 8]
            {
                double y = x * x;
                double z = 2 * Math.Log(x / x3) * bessel_j0(x) / Math.PI;
                r = evaluate_rational(P3, Q3, y);
                factor = (x + x3) * ((x - x31 / 256) - x32);
                value = z + factor * r;
            }
            else                                // x in (8, \infty)
            {
                double y = 8 / x;
                double y2 = y * y;
                double z = x - 0.25f * Math.PI;
                rc = evaluate_rational(PC, QC, y2);
                rs = evaluate_rational(PS, QS, y2);
                factor = Math.Sqrt(2 / (x * Math.PI));
                value = factor * (rc * sin(z) + y * rs * cos(z));
            }

            return value;
        }

        private static double bessel_y1(double x)
        {
            double[] P1 = 
            {
                 4.0535726612579544093e+13,
                 5.4708611716525426053e+12,
                -3.7595974497819597599e+11,
                 7.2144548214502560419e+09,
                -5.9157479997408395984e+07,
                 2.2157953222280260820e+05,
                -3.1714424660046133456e+02,
            };
            double[] Q1 = 
            {
                 3.0737873921079286084e+14,
                 4.1272286200406461981e+12,
                 2.7800352738690585613e+10,
                 1.2250435122182963220e+08,
                 3.8136470753052572164e+05,
                 8.2079908168393867438e+02,
                 1.0,
            };
            double[] P2 = 
            {
                 1.1514276357909013326e+19,
                -5.6808094574724204577e+18,
                -2.3638408497043134724e+16,
                 4.0686275289804744814e+15,
                -5.9530713129741981618e+13,
                 3.7453673962438488783e+11,
                -1.1957961912070617006e+09,
                 1.9153806858264202986e+06,
                -1.2337180442012953128e+03,
            };
            double[] Q2 = 
            {
                 5.3321844313316185697e+20,
                 5.6968198822857178911e+18,
                 3.0837179548112881950e+16,
                 1.1187010065856971027e+14,
                 3.0221766852960403645e+11,
                 6.3550318087088919566e+08,
                 1.0453748201934079734e+06,
                 1.2855164849321609336e+03,
                 1.0,
            };
            double[] PC = 
            {
                -4.4357578167941278571e+06,
                -9.9422465050776411957e+06,
                -6.6033732483649391093e+06,
                -1.5235293511811373833e+06,
                -1.0982405543459346727e+05,
                -1.6116166443246101165e+03,
                 0.0,
            };
            double[] QC = 
            {
                -4.4357578167941278568e+06,
                -9.9341243899345856590e+06,
                -6.5853394797230870728e+06,
                -1.5118095066341608816e+06,
                -1.0726385991103820119e+05,
                -1.4550094401904961825e+03,
                 1.0,
            };
            double[] PS = 
            {
                 3.3220913409857223519e+04,
                 8.5145160675335701966e+04,
                 6.6178836581270835179e+04,
                 1.8494262873223866797e+04,
                 1.7063754290207680021e+03,
                 3.5265133846636032186e+01,
                 0.0,
            };
            double[] QS = 
            {
                 7.0871281941028743574e+05,
                 1.8194580422439972989e+06,
                 1.4194606696037208929e+06,
                 4.0029443582266975117e+05,
                 3.7890229745772202641e+04,
                 8.6383677696049909675e+02,
                 1.0,
            };

            double x1 = 2.1971413260310170351e+00,
                    x2 = 5.4296810407941351328e+00,
                    x11 = 5.620e+02,
                    x12 = 1.8288260310170351490e-03,
                    x21 = 1.3900e+03,
                    x22 = -6.4592058648672279948e-06;

            double value, factor, r, rc, rs;

            if (x <= 0) throw new Exception("x must be > 0");
            if (x <= 4)                       // x in (0, 4]
            {
                double y = x * x;
                double z = 2 * Math.Log(x / x1) * bessel_j1(x) / Math.PI;
                r = evaluate_rational(P1, Q1, y);
                factor = (x + x1) * ((x - x11 / 256) - x12) / x;
                value = z + factor * r;
            }
            else if (x <= 8)                  // x in (4, 8]
            {
                double y = x * x;
                double z = 2 * Math.Log(x / x2) * bessel_j1(x) / Math.PI;
                r = evaluate_rational(P2, Q2, y);
                factor = (x + x2) * ((x - x21 / 256) - x22) / x;
                value = z + factor * r;
            }
            else                                // x in (8, \infty)
            {
                double y = 8 / x;
                double y2 = y * y;
                double z = x - 0.75f * Math.PI;
                rc = evaluate_rational(PC, QC, y2);
                rs = evaluate_rational(PS, QS, y2);
                factor = Math.Sqrt(2 / (x * Math.PI));
                value = factor * (rc * sin(z) + y * rs * cos(z));
            }

            return value;
        }

        #endregion
    }
}
