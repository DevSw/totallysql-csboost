using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class cauchy_distribution : distribution
    {
        double m_a, m_hg;

        public cauchy_distribution(double a, double hg)
        {
            m_a = a;
            m_hg = hg;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (double.IsInfinity(m_a)) throw new ArgumentException(string.Format("Location must be a finite number (got {0:G}).", m_a));
            if (m_hg <= 0 || double.IsInfinity(m_hg)) throw new ArgumentException(string.Format("Scale argument must be a finite number > 0 (got {0:G}).", m_hg));
        }

        public override bool discrete() { return false; }

        public override bool symmetric() { return true; }

        public double location() { return m_a; }

        public double scale() { return m_hg; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(-double.MaxValue, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return range();
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            double xs = (x - m_a) / m_hg;
            return 1 / (Math.PI * m_hg * (1 + xs * xs));
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : -double.MaxValue;
            double x = Math.Sqrt(m_hg / (Math.PI * p) - m_hg * m_hg);
            if (RHS) return m_a + x;
            else return m_a - x;
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            return cdf_imp(x, false);
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            return cdf_imp(x, true);
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            return quantile_imp(p, false);
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            return quantile_imp(q, true);
        }

        public override double mean()
        {
            throw new NotImplementedException("The Cauchy distirbution does not have a mean");
        }

        public override double variance()
        {
            throw new NotImplementedException("The Cauchy distirbution does not have a variance");
        }

        public override double mode()
        {
            return m_a;
        }

        public override double median()
        {
            return m_a;
        }

        public override double skewness()
        {
            throw new NotImplementedException("The Cauchy distirbution does not have a skewness");
        }

        public override double kurtosis_excess()
        {
            throw new NotImplementedException("The Cauchy distirbution does not have a kurtosis");
        }

        public override double standard_deviation()
        {
            throw new NotImplementedException("The Cauchy distirbution does not have a standard deviation");

        }

        private double cdf_imp(double x, bool complement)
        {
            double mx = -Math.Abs((x - m_a) / m_hg); // scale is > 0
            if (mx > -XMath.epsilon / 8) return 0.5;
            double result = -Math.Atan(1 / mx) / Math.PI;
            return (((x > m_a) != complement) ? 1 - result : result);
        }

        private double quantile_imp(double p, bool complement)
        {
            if (p == 1 || p == 0) throw new OverflowException();
            double P = p - Math.Floor(p);   // argument reduction of p:
            if (P > 0.5) P = P - 1;
            if (P == 0.5) return m_a;
            double result = -m_hg / Math.Tan(Math.PI * P);
            return complement ? m_a - result : m_a + result;
        }
    }
    
}
