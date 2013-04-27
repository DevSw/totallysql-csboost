using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class Students_T_Distribution : distribution
    {
        double m_df;

        public Students_T_Distribution(double i)
        {
            m_df = i;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_df <= 0 || double.IsInfinity(m_df)) throw new ArgumentException(string.Format("Degrees of freedom argument must be a finite number > 0 (got {0:G}).", m_df));
        }

        public override bool discrete() { return false; }

        public override bool symmetric() { return true; }

        public double degrees_of_freedom()
        {
            return m_df;
        }

        public override XMath.pair<double> range()
        { // Range of permissible values for random variable x.
            return new XMath.pair<double>(-double.MaxValue, double.MaxValue);
        }

        public override XMath.pair<double> support()
        { // Range of permissible values for random variable x.
            return new XMath.pair<double>(-double.MaxValue, double.MaxValue);
        }

        public override double pdf(double t)
        {
            base.pdf(t);
            double result;
            double basem1 = t * t / m_df;
            if (basem1 < 0.125) result = Math.Exp(-XMath.log1p(basem1) * (1 + m_df) / 2);
            else result = Math.Pow(1 / (1 + basem1), (m_df + 1) / 2);
            result /= Math.Sqrt(m_df) * XMath.beta(m_df / 2, 0.5);
            return result;
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : -double.MaxValue;
            double x = p * Math.Sqrt(m_df) * XMath.beta(m_df / 2, 0.5);
            x = Math.Pow(x, 1.0 / ((m_df + 1) / 2));
            x = Math.Sqrt(m_df / x - m_df);
            if (double.IsInfinity(x)) x = double.MaxValue;
            return RHS ? x : -x;
        }

        public override double cdf(double t)
        {
            base.cdf(t);
            if (t == 0) return 0.5;
            //
            // Calculate probability of Student's t using the incomplete beta function.
            // probability = ibeta(m_df / 2, 1/2, m_df / (m_df + t*t))
            //
            // However when t is small compared to the degrees of freedom, that formula
            // suffers from rounding error, use the identity formula to work around
            // the problem:
            //
            // I[x](a,b) = 1 - I[1-x](b,a)
            //
            // and:
            //
            //     x = df / (df + t^2)
            //
            // so:
            //
            // 1 - x = t^2 / (df + t^2)
            //
            double t2 = t * t;
            double probability;
            if (m_df > 2 * t2)
            {
                double z = t2 / (m_df + t2);
                probability = XMath.ibetac(0.5, m_df / 2, z) / 2;
            }
            else
            {
                double z = m_df / (m_df + t2);
                probability = XMath.ibeta(m_df / 2, 0.5, z) / 2;
            }
            return (t > 0 ? 1 - probability : probability);
        }

        public override double cdfc(double t)
        {
            base.cdfc(t);
            return cdf(-t);
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            double probability = p;

            // Special cases, regardless of m_df.
            if (probability == 0 || probability == 1)
                throw new OverflowException();
            if (probability == 0.5) return 0;
            //
            // Depending on how many digits double has, this may forward
            // to the incomplete beta inverse as above.  Otherwise uses a
            // faster method that is accurate to ~15 digits everywhere
            // and a couple of epsilon at double precision and in the central 
            // region where most use cases will occur...
            //
            return fast_students_t_quantile(m_df, probability);
        }

        public override double quantilec(double p)
        {
            base.quantilec(p);
            return -quantile(p);
        }

        public override double mean()
        {
            return 0;
        }

        public override double variance()
        {
            double v = degrees_of_freedom();
            return v / (v - 2);
        }

        public override double mode()
        {
            return 0;
        }

        public override double median()
        {
            return 0;
        }

        public override double skewness()
        {
            if (m_df <= 3)
                throw new Exception(string.Format("Students-T Distribution: Skewness is undefined for degrees of freedom <= 3, but got {0:G}.", m_df));
            return 0;
        }

        public override double kurtosis()
        {
            if (m_df <= 4)
                throw new Exception(string.Format("Students-T Distribution: Kurtosis is undefined for degrees of freedom <= 4, but got {0:G}.", m_df));
            return 3 * (m_df - 2) / (m_df - 4);
        }

        public override double kurtosis_excess()
        {
            if (m_df <= 4)
                throw new Exception(string.Format("Students-T Distribution: Kurtosis is undefined for degrees of freedom <= 4, but got {0:G}.", m_df));
            return 6 / (m_df - 4);
        }

        private static double slow_students_t_quantile(double df, double p)
        {
            double probability = (p > 0.5) ? 1 - p : p;
            double t, x, y = 0;
            x = XMath.ibeta_inv_imp(df / 2, 0.5, 2 * probability, 1- (2 * probability), ref y);
            if (df * y > double.MaxValue * x) throw new OverflowException();
            else t = Math.Sqrt(df * y / x);
            if (p < 0.5) t = -t;
            return t;
        }

        private static double fast_students_t_quantile(double df, double p)
        {
            bool invert = false;
            if ((df < 2) && (Math.Floor(df) != df)) return slow_students_t_quantile(df, p);
            if (p > 0.5)
            {
                p = 1 - p;
                invert = true;
            }
            //
            // Get an estimate of the result:
            //
            bool exact = false;
            double t = inverse_students_t(df, p, 1 - p, ref exact);
            if ((t == 0) || exact)
                return invert ? -t : t; // can't do better!
            //
            // Change variables to inverse incomplete beta:
            //
            double t2 = t * t;
            double xb = df / (df + t2);
            double y = t2 / (df + t2);
            double a = df / 2;
            //
            // t can be so large that x underflows,
            // just return our estimate in that case:
            //
            if (xb == 0)
                return t;
            //
            // Get incomplete beta and it's derivative:
            //
            double f1 = 0;
            double f0 = xb < y ? XMath.ibeta_imp(a, 0.5, xb, false, true, ref f1) : XMath.ibeta_imp(0.5, a, y, true, true, ref f1);

            // Get cdf from incomplete beta result:
            double p0 = f0 / 2 - p;
            // Get pdf from derivative:
            double p1 = f1 * Math.Sqrt(y * xb * xb * xb / df);
            //
            // Second derivative divided by p1:
            //
            // yacas gives:
            //
            // In> PrettyForm(Simplify(D(t) (1 + t^2/v) ^ (-(v+1)/2)))
            //
            //  |                        | v + 1     |     |
            //  |                       -| ----- + 1 |     |
            //  |                        |   2       |     |
            // -|             |  2     |                   |
            //  |             | t      |                   |
            //  |             | -- + 1 |                   |
            //  | ( v + 1 ) * | v      |               * t |
            // ---------------------------------------------
            //                       v
            //
            // Which after some manipulation is:
            //
            // -p1 * t * (df + 1) / (t^2 + df)
            //
            double p2 = t * (df + 1) / (t * t + df);
            // Halley step:
            t = Math.Abs(t);
            t += p0 / (p1 + p0 * p2 / 2);
            return !invert ? -t : t;
        }

        public static double inverse_students_t(double df, double u, double v, ref bool pexact)
        {
            //
            // df = number of degrees of freedom.
            // u = probablity.
            // v = 1 - u.
            //
            bool invert = false;
            double result = 0;
            pexact = false;
            if (u > v)
            {
                // function is symmetric, invert it:
                XMath.swap(ref u, ref v);
                invert = true;
            }
            if (df == 1.0 || df == 2.0 || df == 4.0 || df == 6.0)
            {
                //
                // we have integer degrees of freedom, try for the special
                // cases first:
                //
                switch ((int)df)
                {
                    case 1:
                        {
                            //
                            // df = 1 is the same as the Cauchy distribution, see
                            // Shaw Eq 35:
                            //
                            if (u == 0.5)
                                result = 0;
                            else
                                result = -Math.Cos(Math.PI * u) / Math.Sin(Math.PI * u);
                            pexact = true;
                            break;
                        }
                    case 2:
                        {
                            //
                            // df = 2 has an exact result, see Shaw Eq 36:
                            //
                            result = (2 * u - 1) / Math.Sqrt(2 * u * v);
                            pexact = true;
                            break;
                        }
                    case 4:
                        {
                            //
                            // df = 4 has an exact result, see Shaw Eq 38 & 39:
                            //
                            double alpha = 4 * u * v;
                            double root_alpha = Math.Sqrt(alpha);
                            double r = 4 * Math.Cos(Math.Acos(root_alpha) / 3) / root_alpha;
                            double x = Math.Sqrt(r - 4);
                            result = u - 0.5f < 0 ? (double)-x : x;
                            pexact = true;
                            break;
                        }
                    case 6:
                        {
                            double tolerance = XMath.ldexp(1.0, (2 * 53) / 3);
                            //
                            // We get numeric overflow in this area:
                            //
                            if (u < 1e-150)
                                return (invert ? -1 : 1) * inverse_students_t_hill(df, u);
                            //
                            // Newton-Raphson iteration of a polynomial case,
                            // choice of seed value is taken from Shaw's online
                            // supplement:
                            //
                            double a = 4 * (u - u * u);//1 - 4 * (u - 0.5f) * (u - 0.5f);
                            double b = Math.Pow(a, 1.0 / 3.0);                                  //CBRT
                            double c = 0.85498797333834849467655443627193;
                            double p = 6 * (1 + c * (1 / b - 1));
                            double p0;
                            do
                            {
                                double p2 = p * p;
                                double p4 = p2 * p2;
                                double p5 = p * p4;
                                p0 = p;
                                // next term is given by Eq 41:
                                p = 2 * (8 * a * p5 - 270 * p2 + 2187) / (5 * (4 * a * p4 - 216 * p - 243));
                            } while (Math.Abs((p - p0) / p) > tolerance);
                            //
                            // Use Eq 45 to extract the result:
                            //
                            p = Math.Sqrt(p - df);
                            result = (u - 0.5) < 0 ? -p : p;
                            break;
                        }
                }
            }
            else
            {
                if (df < 3)
                {
                    //
                    // Use a roughly linear scheme to choose between Shaw's
                    // tail series and body series:
                    //
                    double crossover = 0.2742 - df * 0.0242143;
                    if (u > crossover)
                    {
                        result = inverse_students_t_body_series(df, u);
                    }
                    else
                    {
                        result = inverse_students_t_tail_series(df, u);
                    }
                }
                else
                {
                    //
                    // Use Hill's method except in the exteme tails
                    // where we use Shaw's tail series.
                    // doublehe crossover point is roughly exponential in -df:
                    //
                    double crossover = XMath.ldexp(1.0, (int)Math.Round(df / -0.654));
                    if (u > crossover)
                    {
                        result = inverse_students_t_hill(df, u);
                    }
                    else
                    {
                        result = inverse_students_t_tail_series(df, u);
                    }
                }
            }
            return invert ? (double)-result : result;
        }

        public static double find_ibeta_inv_from_t_dist(double a, double p, double q, ref double py)
        {
            double u = (p > q) ? (0.5 - q) / 2.0 : p / 2.0;
            double v = 1.0 - u; // u < 0.5 so no cancellation error
            double df = a * 2.0;
            bool pexact = false;
            double t = inverse_students_t(df, u, v, ref pexact);
            double x = df / (df + t * t);
            py = t * t / (df + t * t);
            return x;
        }

        private static double inverse_students_t_hill(double ndf, double u)
        {
            double a, b, c, d, q, x, y;

            double root_two = XMath.root_two;
            if (ndf > 1e20) return -XMath.erfc_inv(2 * u) * root_two;

            a = 1 / (ndf - 0.5);
            b = 48 / (a * a);
            c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
            d = ((94.5 / (b + c) - 3) / b + 1) * Math.Sqrt(a * Math.PI / 2.0) * ndf;
            y = Math.Pow(d * 2 * u, 2 / ndf);

            if (y > (0.05 + a))
            {
                //
                // Asymptotic inverse expansion about normal:
                //
                x = -XMath.erfc_inv(2 * u) * root_two;
                y = x * x;

                if (ndf < 5)
                    c += 0.3 * (ndf - 4.5) * (x + 0.6);
                c += (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b;
                y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x;
                y = XMath.expm1(a * y * y);
            }
            else
            {
                y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822)
                        * (ndf + 2) * 3) + 0.5 / (ndf + 4)) * y - 1)
                        * (ndf + 1) / (ndf + 2) + 1 / y;
            }
            q = Math.Sqrt(ndf * y);

            return -q;
        }

        private static double inverse_students_t_tail_series(double df, double v)
        {
            // Tail series expansion, see section 6 of Shaw's paper.
            // w is calculated using Eq 60:
            double w = XMath.gamma_delta_ratio(df / 2, 0.5)
               * Math.Sqrt(df * Math.PI) * v;
            // define some variables:
            double np2 = df + 2;
            double np4 = df + 4;
            double np6 = df + 6;
            //
            // Calculate the coefficients d(k), these depend only on the
            // number of degrees of freedom df, so at least in theory
            // we could tabulate these for fixed df, see p15 of Shaw:
            //
            double[] d = new double[7];
            d[0] = 1.0;
            d[1] = -(df + 1) / (2 * np2);
            np2 *= (df + 2);
            d[2] = -df * (df + 1) * (df + 3) / (8 * np2 * np4);
            np2 *= df + 2;
            d[3] = -df * (df + 1) * (df + 5) * (((3 * df) + 7) * df - 2) / (48 * np2 * np4 * np6);
            np2 *= (df + 2);
            np4 *= (df + 4);
            d[4] = -df * (df + 1) * (df + 7) *
               ((((((15 * df) + 154) * df + 465) * df + 286) * df - 336) * df + 64)
               / (384 * np2 * np4 * np6 * (df + 8));
            np2 *= (df + 2);
            d[5] = -df * (df + 1) * (df + 3) * (df + 9)
                     * (((((((35 * df + 452) * df + 1573) * df + 600) * df - 2020) * df) + 928) * df - 128)
                     / (1280 * np2 * np4 * np6 * (df + 8) * (df + 10));
            np2 *= (df + 2);
            np4 *= (df + 4);
            np6 *= (df + 6);
            d[6] = -df * (df + 1) * (df + 11)
                     * ((((((((((((945 * df) + 31506) * df + 425858) * df + 2980236) * df + 11266745) * df + 20675018) * df + 7747124) * df - 22574632) * df - 8565600) * df + 18108416) * df - 7099392) * df + 884736)
                     / (46080 * np2 * np4 * np6 * (df + 8) * (df + 10) * (df + 12));
            //
            // Now bring everthing together to provide the result,
            // this is Eq 62 of Shaw:
            //
            double rn = Math.Sqrt(df);
            double div = Math.Pow(rn * w, 1 / df);
            double power = div * div;
            double result = XMath.evaluate_polynomial(d, power);
            result *= rn;
            result /= div;
            return -result;
        }

        private static double inverse_students_t_body_series(double df, double u)
        {
            //
            // Body series for small N:
            // Start with Eq 56 of Shaw:
            //
            double v = XMath.gamma_delta_ratio(df / 2, 0.5)
               * Math.Sqrt(df * Math.PI) * (u - 0.5);
            //
            // Workspace for the polynomial coefficients:
            //
            double[] c = new double[11];
            c[0] = 0.0;
            c[1] = 1.0;
            //
            // Figure out what the coefficients are, note these depend
            // only on the degrees of freedom (Eq 57 of Shaw):
            //
            double cf = 1 / df;
            c[2] = 0.16666666666666666667 + 0.16666666666666666667 * cf;
            c[3] = (0.0083333333333333333333 * cf
               + 0.066666666666666666667) * cf
               + 0.058333333333333333333;
            c[4] = ((0.00019841269841269841270 * cf
               + 0.0017857142857142857143) * cf
               + 0.026785714285714285714) * cf
               + 0.025198412698412698413;
            c[5] = (((2.7557319223985890653e-6 * cf
               + 0.00037477954144620811287) * cf
               - 0.0011078042328042328042) * cf
               + 0.010559964726631393298) * cf
               + 0.012039792768959435626;
            c[6] = ((((2.5052108385441718775e-8 * cf
               - 0.000062705427288760622094) * cf
               + 0.00059458674042007375341) * cf
               - 0.0016095979637646304313) * cf
               + 0.0061039211560044893378) * cf
               + 0.0038370059724226390893;
            c[7] = (((((1.6059043836821614599e-10 * cf
               + 0.000015401265401265401265) * cf
               - 0.00016376804137220803887) * cf
               + 0.00069084207973096861986) * cf
               - 0.0012579159844784844785) * cf
               + 0.0010898206731540064873) * cf
               + 0.0032177478835464946576;
            c[8] = ((((((7.6471637318198164759e-13 * cf
               - 3.9851014346715404916e-6) * cf
               + 0.000049255746366361445727) * cf
               - 0.00024947258047043099953) * cf
               + 0.00064513046951456342991) * cf
               - 0.00076245135440323932387) * cf
               + 0.000033530976880017885309) * cf
               + 0.0017438262298340009980;
            c[9] = (((((((2.8114572543455207632e-15 * cf
               + 1.0914179173496789432e-6) * cf
               - 0.000015303004486655377567) * cf
               + 0.000090867107935219902229) * cf
               - 0.00029133414466938067350) * cf
               + 0.00051406605788341121363) * cf
               - 0.00036307660358786885787) * cf
               - 0.00031101086326318780412) * cf
               + 0.00096472747321388644237;
            c[10] = ((((((((8.2206352466243297170e-18 * cf
               - 3.1239569599829868045e-7) * cf
               + 4.8903045291975346210e-6) * cf
               - 0.000033202652391372058698) * cf
               + 0.00012645437628698076975) * cf
               - 0.00028690924218514613987) * cf
               + 0.00035764655430568632777) * cf
               - 0.00010230378073700412687) * cf
               - 0.00036942667800009661203) * cf
               + 0.00054229262813129686486;
            //
            // The result is then a polynomial cf v (see Eq 56 of Shaw):
            //
            return XMath.evaluate_odd_polynomial(c, v);
        }

    }
}
