using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class normal_distribution : distribution
    {
        double m_mean;  // distribution mean or location.
        double m_sd;    // distribution standard deviation or scale.
        double root_two = XMath.root_two;

        public normal_distribution(double mean, double sd)
        {
            m_mean = mean;
            m_sd = sd;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (double.IsInfinity(m_mean)) throw new ArgumentException(string.Format("Mean must be a finite number (got {0:G}).", m_mean));
            if (m_sd <= 0 || double.IsInfinity(m_sd)) throw new ArgumentException(string.Format("Standard deviation must be a finite number > 0 (got {0:G}).", m_sd));
        }

        public override bool discrete() { return false; }

        public override bool symmetric() { return true; }

        public override double mean()
        { // alias for location.
            return m_mean;
        }

        public override double variance()
        {
            return m_sd * m_sd;
        }

        public override double standard_deviation()
        { // alias for scale.
            return m_sd;
        }

        // Synonyms, provided to allow generic use of find_location and find_scale.
        public double location()
        { // location.
            return m_mean;
        }

        public double scale()
        { // scale.
            return m_sd;
        }

        public override XMath.pair<double> range()
        { // Range of permissible values for random variable x.
            return new XMath.pair<double>(-double.MaxValue, double.MaxValue); // - to + max value.
        }

        public override XMath.pair<double> support()
        { // Range of supported values for random variable x.
            // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
            return new XMath.pair<double>(-double.MaxValue, double.MaxValue); // - to + max value.
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            if (double.IsInfinity(x)) return 0; // pdf + and - infinity is zero.

            double exponent = x - m_mean;
            exponent *= -exponent;
            exponent /= 2 * m_sd * m_sd;

            double result = Math.Exp(exponent);
            result /= m_sd * XMath.root_two_pi;

            return result;
        } // pdf

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : -double.MaxValue;
            double x = Math.Sqrt(-Math.Log(p * m_sd * XMath.root_two_pi) * 2 * m_sd * m_sd);
            if (RHS) return m_mean + x;
            return m_mean - x;
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            double result;
            if (double.IsInfinity(x))
            {
                if (x < 0) return 0; // -infinity
                return 1; // + infinity
            }
            double diff = (x - m_mean) / (m_sd * root_two);
            result = XMath.erfc(-diff) / 2.0;
            return result;
        } // cdf

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (double.IsInfinity(x))
            {
                if (x < 0) return 1; // cdf complement -infinity is unity.
                return 0; // cdf complement +infinity is zero
            }
            double result;
            double diff = (x - m_mean) / (m_sd * root_two);
            result = XMath.erfc(diff) / 2.0;
            return result;
        } // cdf complement

        public override double quantile(double p)
        {
            base.quantile(p);
            double result;
            result = XMath.erfc_inv(2 * p);
            result = -result;
            result *= m_sd * root_two;
            result += m_mean;
            return result;
        } // quantile

        public override double quantilec(double q)
        {
            base.quantilec(q);
            double result;
            result = XMath.erfc_inv(2 * q);
            result *= m_sd * root_two;
            result += m_mean;
            return result;
        } // quantile

        public override double mode()
        {
            return mean();
        }

        public override double median()
        {
            return mean();
        }

        public override double skewness()
        {
            return 0;
        }

        public override double kurtosis()
        {
            return 3;
        }

        public override double kurtosis_excess()
        {
            return 0;
        }

    } // class normal_distribution
}
