using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class inverse_chisquared_distribution : distribution
    {
        double m_df, m_scale;

        public inverse_chisquared_distribution(double df, double scale)
        {
            m_df = df;
            m_scale = scale;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_df <= 0 || double.IsInfinity(m_df)) throw new ArgumentException(string.Format("Inverse Chi-Squared Distribution: Degrees of freedom argument must be a finite number > 0 (got {0:G}).", m_df));
            if (m_scale <= 0 || double.IsInfinity(m_scale)) throw new ArgumentException(string.Format("Inverse Chi-Squared Distribution: Scale argument must be a finite number > 0 (got {0:G}).", m_scale));
        }

        public override bool discrete() { return false; }

        public double degrees_of_freedom() { return m_df; }

        public double scale() { return m_scale; }

        public override bool LHS() { return true; }

        public override bool tail_left() { return false; }

        public override bool unimodal()
        {
            return true;
        }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            if (x == 0) return 0;

            double result = (m_df * m_scale/2) / x;
            if(result < double.Epsilon) return 0; // Random variable is near enough infinite.
            result = XMath.gamma_p_derivative(m_df/2, result) * m_df * m_scale/2;
            if(result != 0) result /= (x * x);
            return result;
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : 0;
            if (p == max_pdf()) return mode();

            double ubound = RHS ? double.MaxValue : mode();
            double lbound = RHS ? mode() : 0;
            return find_pdf_inv(p, lbound, ubound, RHS ? false : true);
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if (x == 0) return 0;
            return XMath.gamma_q(m_df / 2, (m_df * (m_scale / 2)) / x);        
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (x == 0) return 1;
            return XMath.gamma_p(m_df / 2, (m_df * m_scale / 2) / x); 
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            double result = XMath.gamma_q_inv(m_df /2, p);
            if(result == 0) return double.PositiveInfinity;
            result = m_df * (m_scale / 2) / result;
            return result;
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            double result = XMath.gamma_p_inv(m_df / 2, q);
            if (result == 0) return double.PositiveInfinity;
            result = (m_df * m_scale / 2) / result;
            return result;
        }

        public override double mean()
        {
            if (m_df <= 2) throw new Exception(string.Format("Inverse Chi Squared Distribution: The mean is only defined when df > 2 (got {0:G}).", m_df));
            return (m_df * m_scale) / (m_df - 2);
        }

        public override double variance()
        {
            if (m_df <= 4) throw new Exception(string.Format("Inverse Chi Squared Distribution: The variance is only defined when df > 4 (got {0:G}).", m_df));
            return 2 * m_df * m_df * m_scale * m_scale / ((m_df - 2) * (m_df - 2) * (m_df - 4));
        }

        public override double mode()
        {
            return (m_df * m_scale) / (m_df + 2);
        }

        //Median supplied by base class

        public override double skewness()
        {
            if (m_df <= 6) throw new Exception(string.Format("Inverse Chi Squared Distribution: The skewness is only defined when df > 6 (got {0:G}).", m_df));
            return 4 * Math.Sqrt(2 * (m_df - 4)) / (m_df - 6);  // Not a function of scale.
        }

        //Kurtosis supplied by base class

        public override double kurtosis_excess()
        {
            if (m_df <= 8) throw new Exception(string.Format("Inverse Chi Squared Distribution: The kurtosis is only defined when df > 8 (got {0:G}).", m_df));
            return 12 * (5 * m_df - 22) / ((m_df - 6) * (m_df - 8));  // Not a function of scale.
        }
    }
}
