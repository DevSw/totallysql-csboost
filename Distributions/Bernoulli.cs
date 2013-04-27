using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class bernoulli_distribution: distribution
    {
        double m_p;

        public bernoulli_distribution(double p)
        {
            m_p = p;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_p < 0 || m_p > 1) throw new ArgumentException(string.Format("Success fraction must be >= 0 and <=1 (got {0:G}).", m_p));        
        }

        public override bool discrete() { return true; }

        public double success_fraction() { return m_p; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, 1);
        }

        public override XMath.pair<double> support()
        {
            return range();
        }

        public override void setup_bars()
        {
            return;
        }

        public override bool LHS() { return m_p > 0.5; }

        public override bool RHS() { return m_p < 0.5; }

        public override bool tail_right() { return false; }

        public override bool tail_left() { return false; }

        public override double mean()
        {
            return m_p;
        }

        public override double variance()
        {
            return m_p * (1 - m_p);
        }

        public override double max_pdf()
        {
            if (m_p < 0.5) return pdf(0);
            return pdf(1);
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            if (x == 0) return 1 - m_p;
            return m_p;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            return p == m_p ? 1 : 0;
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if (x == 0) return 1 - m_p;
            return 1;
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (x == 0) return m_p;
            return 0;
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p <= 1 - m_p) return 0;
            return 1;
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            if (q <= 1 - m_p) return 1;
            return 0;
        }

        public override double mode()
        {
            return m_p <= 0.5 ? 0 : 1;
        }

        public override double skewness()
        {
            return (1 - 2 * m_p) / Math.Sqrt(m_p * (1 - m_p));
        }

        public override double kurtosis_excess()
        {
            return 1 / (1 - m_p) + 1 / m_p - 6;
        }
        public override double kurtosis()
        {
            return 1 / (1 - m_p) + 1 / m_p - 6 + 3;
        }

   }
}
