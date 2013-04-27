using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class chisquared_distribution : distribution
    {
        double m_df;

        public chisquared_distribution(double i)
        {
            m_df = i;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_df <= 0 || double.IsInfinity(m_df)) throw new ArgumentException(string.Format("Degrees of freedom argument must be a finite number > 0 (got {0:G}).", m_df));
        }

        public override bool discrete() { return false; }

        public double degrees_of_freedom() { return m_df; }

        public override bool LHS() { return m_df > 2; }

        public override bool DCL()
        {
            return m_df < 2;
        }

        public override bool tail_left() { return false; }

        public override bool unimodal()
        {
            return m_df > 2;
        }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(m_df > +2 ? 0 : 2 * double.Epsilon, double.MaxValue);
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            if (x == 0)
            {
                // Handle special cases:
                if (m_df < 2) return double.MaxValue;
                else if (m_df == 2) return 0.5;
                else return 0;
            }
            return XMath.gamma_p_derivative(m_df / 2, x / 2) / 2.0;
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double max_pdf()
        {
            if (m_df == 2) return pdf(0);
            if (m_df < 2) return double.MaxValue;
            return pdf(mode());
        }

        public override double pdf_inv(double p, bool RHS)
        {
            if (m_df == 2 && !RHS) throw new Exception("The Chi Squared Distribution has no LHS when degrees of freedom = 2 - so no LHS inverse pdf is available.");
            base.pdf_inv(p, RHS);
            if (p == 0) return m_df == 2 || RHS ? double.MaxValue : 0;
            double ubound = quantilec(2 * double.Epsilon);
            double lbound = DCL() ? quantile(XMath.epsilon / m_df) : 0;
            if (p < 1e-300) ubound *= 1.5;                                                 //Workaround for unpredictable results when p is very small
            if (RHS) return find_pdf_inv(p, unimodal() ? mode() : lbound, ubound, false);
            else return find_pdf_inv(p, lbound, unimodal() ? mode() : ubound, true);
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            return XMath.gamma_p(m_df / 2, x / 2);
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            return XMath.gamma_q(m_df / 2, x / 2);
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            return 2.0 * XMath.gamma_p_inv(m_df / 2, p);
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            return 2.0 * XMath.gamma_q_inv(m_df / 2, q);
        }

        public override double mean()
        {
            return m_df;
        }

        public override double variance()
        {
            return 2.0 * m_df;
        }

        public override double mode()
        {
            if (m_df < 2) throw new Exception(string.Format("ChiSquare.mode: mode is only defined when degrees of freedom are >= 2 (DF = {0:G}).", m_df));
            return m_df - 2;
        }

        //Median supplied by base class

        public override double skewness()
        {
            return Math.Sqrt(8.0 / m_df);
        }

        //Kurtosis supplied by base class

        public override double kurtosis_excess()
        {
            return 12.0 / m_df;
        }
    }
}
