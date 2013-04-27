using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{

    public class fisher_f_distribution : distribution
    {
        double m_df1;
        double m_df2;

        public fisher_f_distribution(double i, double j)
        {
            m_df1 = i;
            m_df2 = j;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_df1 <= 0 || double.IsInfinity(m_df1)) throw new ArgumentException(string.Format("Degrees of freedom argument must be a finite number > 0 (got {0:G}).", m_df1));
            if (m_df2 <= 0 || double.IsInfinity(m_df2)) throw new ArgumentException(string.Format("Degrees of freedom argument must be a finite number > 0 (got {0:G}).", m_df2));
        }

        public override bool discrete() { return false; }

        double degrees_of_freedom1()
        {
            return m_df1;
        }

        double degrees_of_freedom2()
        {
            return m_df2;
        }

        public override bool LHS()
        {
            return m_df1 > 2;
        }

        public override bool DCL()
        {
            return m_df1 < 2;
        }

        public override bool tail_left() { return false; }

        public override bool unimodal()
        {
            return m_df1 > 2;
        }

        public override XMath.pair<double> range()
        { // Range of permissible values for random variable x.
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        { // Range of permissible values for random variable x.
            return new XMath.pair<double>(DCL() ? double.Epsilon: 0, double.MaxValue);
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();

            if ((x < 0) || double.IsPositiveInfinity(x))
                throw new Exception(string.Format("Random variable parameter was {0:G}, but must be > 0 !", x));

            if (x == 0)
            {
                // special cases:
                if (df1 < 2)
                    throw new OverflowException();
                else if (df1 == 2)
                    return 1;
                else
                    return 0;
            }

            //
            // You reach this formula by direct differentiation of the
            // cdf expressed in terms of the incomplete beta.
            //
            // There are two versions so we don't pass a value of z
            // that is very close to 1 to ibeta_derivative: for some values
            // of df1 and df2, all the change takes place in this area.
            //
            double v1x = df1 * x;
            if (x > 0 && df1 > 0 && v1x < 2 * double.Epsilon) return double.MaxValue;   //avoid the overflow error
            if (double.IsPositiveInfinity(v1x)) return 0;
            double result;
            if (v1x > df2)
            {
                result = (df2 * df1) / ((df2 + v1x) * (df2 + v1x));
                result *= XMath.ibeta_derivative(df2 / 2, df1 / 2, df2 / (df2 + v1x));
            }
            else
            {
                result = df2 + df1 * x;
                result = (result * df1 - x * df1 * df1) / (result * result);
                result *= XMath.ibeta_derivative(df1 / 2, df2 / 2, v1x / (df2 + v1x));
            }
            return result;
        } // pdf

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            double lbound, ubound;
            if (p == 0) return RHS ? double.MaxValue : 0;
            if (RHS)
            {
                lbound = unimodal() ? mode() : support().v1;
                ubound = !unimodal() && p > 1 ? 1.0 / p : support().v2;
            }
            else
            {
                lbound = support().v1;
                ubound = mode();            //must be unimodal or we wouldn't be here
            }
            return find_pdf_inv(p, lbound, ubound, !RHS);     
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();

            if ((x < 0) || double.IsPositiveInfinity(x))
                throw new Exception(string.Format("Random variable parameter was {0:G}, but must be > 0 !", x));

            double v1x = df1 * x;
            //
            // There are two equivalent formulas used here, the aim is
            // to prevent the final argument to the incomplete beta
            // from being too close to 1: for some values of df1 and df2
            // the rate of change can be arbitrarily large in this area,
            // whilst the value we're passing will have lost information
            // content as a result of being 0.999999something.  Better
            // to switch things around so we're passing 1-z instead.
            //
            return v1x > df2 ? XMath.ibetac(df2 / 2, df1 / 2, df2 / (df2 + v1x)) : XMath.ibeta(df1 / 2, df2 / 2, v1x / (df2 + v1x));
        } // cdf

        public override double cdfc(double x)
        {
            base.cdfc(x);
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();

            if ((x < 0) || double.IsPositiveInfinity(x))
                throw new Exception(string.Format("Random variable parameter was {0:G}, but must be > 0 !", x));

            double v1x = df1 * x;
            //
            // There are two equivalent formulas used here, the aim is
            // to prevent the final argument to the incomplete beta
            // from being too close to 1: for some values of df1 and df2
            // the rate of change can be arbitrarily large in this area,
            // whilst the value we're passing will have lost information
            // content as a result of being 0.999999something.  Better
            // to switch things around so we're passing 1-z instead.
            //
            return v1x > df2 ? XMath.ibeta(df2 / 2, df1 / 2, df2 / (df2 + v1x)) : XMath.ibetac(df1 / 2, df2 / 2, v1x / (df2 + v1x));
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();

            double x, y = 0.0;

            x = XMath.ibeta_inv_imp(df1 / 2, df2 / 2, p, 1.0 - p, ref y);

            double result = df2 * x / (df1 * y);
            if (result == 0 && p > 0) result = double.Epsilon;
            return result;
        } // quantile

        public override double quantilec(double q)
        {
            base.quantilec(q);
            double x, y = 0;
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();

            x = XMath.ibeta_inv_imp(df1 / 2, df2 / 2, 1-q, q, ref y);

            double result = df2 * x / (df1 * y);
            return result > double.MaxValue ? double.MaxValue : result;
        }

        public override double mean()
        { // Mean of F distribution = v.
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();

            if (df2 <= 2)
                throw new Exception(string.Format("Second degree of freedom was {0:G} but must be > 2 in order for the distribution to have a mean.", df2));
            return df2 / (df2 - 2);
        } // mean

        public override double variance()
        { // Variance of F distribution.
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();

            if (df2 <= 4)
                throw new Exception(string.Format("Second degree of freedom was {0:G} but must be > 4 in order for the distribution to have a valid variance.", df2));
            return 2 * df2 * df2 * (df1 + df2 - 2) / (df1 * (df2 - 2) * (df2 - 2) * (df2 - 4));
        } // variance

        public override double mode()
        {
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();
            if(df1 <= 2) throw new Exception(string.Format("F Distribution: mode is only defined when the first degree of freedom > 2 (got {0:G}.", df1));
            return df2 * (df1 - 2) / (df1 * (df2 + 2));
        }

        //Median supplied by base class

        public override double skewness()
        {
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();

            if (df2 <= 6)
                throw new Exception(string.Format("Second degree of freedom was {0:G} but must be > 6 in order for the distribution to have a skewness.", df2));
            return 2 * (df2 + 2 * df1 - 2) * Math.Sqrt((2 * df2 - 8) / (df1 * (df2 + df1 - 2))) / (df2 - 6);
        }

        //Kurtosis supplied by basis class

        public override double kurtosis_excess()
        {
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();
            if (df2 <= 8)
                throw new Exception(string.Format("Second degree of freedom was {0:G} but must be > 8 in order for the distribution to have a kurtosis.", df2));
            double df2_2 = df2 * df2;
            double df1_2 = df1 * df1;
            double n = -16 + 20 * df2 - 8 * df2_2 + df2_2 * df2 + 44 * df1 - 32 * df2 * df1 + 5 * df2_2 * df1 - 22 * df1_2 + 5 * df2 * df1_2;
            n *= 12;
            double d = df1 * (df2 - 6) * (df2 - 8) * (df1 + df2 - 2);
            return n / d;
        }
    }
}

