using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{

    public class non_central_f_distribution : distribution
    {
        double m_df1;
        double m_df2;
        double m_lambda;

        public non_central_f_distribution(double i, double j, double lambda)
        {
            m_df1 = i;
            m_df2 = j;
            m_lambda = lambda;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_df1 <= 0 || double.IsInfinity(m_df1)) throw new ArgumentException(string.Format("Degrees of freedom argument must be a finite number > 0 (got {0:G}).", m_df1));
            if (m_df2 <= 0 || double.IsInfinity(m_df2)) throw new ArgumentException(string.Format("Degrees of freedom argument must be a finite number > 0 (got {0:G}).", m_df2));
            if (m_lambda < 0 || double.IsInfinity(m_lambda)) throw new ArgumentException(string.Format("Non-centrality argument must be a finite number >= 0 (got {0:G}).", m_lambda));
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

        double non_centrality()
        {
            return m_lambda;
        }

        public override bool LHS()
        {
            return m_df1 > 2;
        }

        public override bool DCL()
        {
            return false;
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
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            double a = m_df1 / 2;
            double b = m_df2 / 2;
            double y = x * a / b;
            double result = new non_central_beta_distribution(a, b, m_lambda).pdf(y / (1 + y));
            result *= (m_df1 / m_df2) / ((1 + y) * (1 + y));
            return result;
        }

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

         double alpha = degrees_of_freedom1() / 2;
         double beta = degrees_of_freedom2() / 2;
         double y = x * alpha / beta;
         double c = y / (1 + y);
         double cp = 1 / (1 + y);
         //
         // To ensure accuracy, we pass both x and 1-x to the
         // non-central beta cdf routine, this ensures accuracy
         // even when we compute x to be ~ 1:
         //
         double r = new non_central_beta_distribution(alpha, beta, m_lambda).non_central_beta_cdf(c, cp, alpha, beta, m_lambda, false);
         return r;
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();

            if ((x < 0) || double.IsPositiveInfinity(x))
                throw new Exception(string.Format("Random variable parameter was {0:G}, but must be > 0 !", x));

            double alpha = degrees_of_freedom1() / 2;
            double beta = degrees_of_freedom2() / 2;
            double y = x * alpha / beta;
            double c = y / (1 + y);
            double cp = 1 / (1 + y);
            //
            // To ensure accuracy, we pass both x and 1-x to the
            // non-central beta cdf routine, this ensures accuracy
            // even when we compute x to be ~ 1:
            //
            double r = new non_central_beta_distribution(alpha, beta, m_lambda).non_central_beta_cdf(c, cp, alpha, beta, m_lambda, true);
            return r;
        }

        public override double quantile(double p)
        {
            base.quantile(p);

            double alpha = m_df1 / 2;
            double beta = m_df2 / 2;
            double x = new non_central_beta_distribution(alpha, beta, m_lambda).quantile(p);

            if (x == 1) throw new OverflowException();
            return (x / (1 - x)) * m_df2 / m_df1;
        } // quantile

        public override double quantilec(double q)
        {
            base.quantilec(q);

            double alpha = m_df1 / 2;
            double beta = m_df2 / 2;
            double x = new non_central_beta_distribution(alpha, beta, m_lambda).quantilec(q);

            if (x == 1) throw new OverflowException();
            return (x / (1 - x)) * (m_df2 / m_df1);
        }

        public override double mean()
        { // Mean of F distribution = v.

            if (m_df2 <= 2)
                throw new Exception(string.Format("Second degree of freedom was {0:G} but must be > 2 in order for the distribution to have a mean.", m_df2));
            return m_df2 * (m_df1 + m_lambda) / (m_df1 * (m_df2 - 2));

        } // mean

        public override double variance()
        { // Variance of F distribution.
            double n = degrees_of_freedom1();
            double m = degrees_of_freedom2();
            double l = non_centrality();

            if (m <= 4)
                throw new Exception(string.Format("Second degree of freedom was {0:G} but must be > 4 in order for the distribution to have a valid variance.", m));
            double result = 2 * m * m * ((n + l) * (n + l)
               + (m - 2) * (n + 2 * l));
            result /= (m - 4) * (m - 2) * (m - 2) * n * n;
            return result;
        } // variance

        public override double mode()
        {
            double df1 = degrees_of_freedom1();
            double df2 = degrees_of_freedom2();
            if(df1 <= 2) throw new Exception(string.Format("Non-Central F Distribution: mode is only defined when the first degree of freedom > 2 (got {0:G}.", df1));
            return generic_find_mode(df2 * (df1 + m_lambda) / (df1 * (df2 - 2)), 0);
        }

        //Median supplied by base class

        public override double skewness()
        {
            double n = degrees_of_freedom1();
            double m = degrees_of_freedom2();
            double l = non_centrality();

            if (m <= 6)
                throw new Exception(string.Format("Second degree of freedom was {0:G} but must be > 6 in order for the distribution to have a skewness.", m));
         double result = 2 * XMath.root_two;
         result *= Math.Sqrt(m - 4);
         result *= (n * (m + n - 2) *(m + 2 * n - 2)
            + 3 * (m + n - 2) * (m + 2 * n - 2) * l
            + 6 * (m + n - 2) * l * l + 2 * l * l * l);
         result /= (m - 6) * Math.Pow(n * (m + n - 2) + 2 * (m + n - 2) * l + l * l, 1.5);
         return result;
        }

        //Kurtosis supplied by basis class

        public override double kurtosis_excess()
        {
            double n = degrees_of_freedom1();
            double m = degrees_of_freedom2();
            double l = non_centrality();

            if (m <= 8)
                throw new Exception(string.Format("Second degree of freedom was {0:G} but must be > 8 in order for the distribution to have a kurtosis.", m));
            double l2 = l * l;
            double l3 = l2 * l;
            double l4 = l2 * l2;
            double result = (3 * (m - 4) * (n * (m + n - 2)
                * (4 * (m - 2) * (m - 2)
                + (m - 2) * (m + 10) * n
                + (10 + m) * n * n)
                + 4 * (m + n - 2) * (4 * (m - 2) * (m - 2)
                + (m - 2) * (10 + m) * n
                + (10 + m) * n * n) * l + 2 * (10 + m)
                * (m + n - 2) * (2 * m + 3 * n - 4) * l2
                + 4 * (10 + m) * (-2 + m + n) * l3
                + (10 + m) * l4))
                /
                ((-8 + m) * (-6 + m) * Math.Pow(n * (-2 + m + n) + 2 * (-2 + m + n) * l + l2, 2));
            return result;
        }
    }
}

