using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class rayleigh_distribution : distribution
    {
        double m_sigma;

        public rayleigh_distribution(double sigma)
        {
            m_sigma = sigma;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_sigma <= 0 || double.IsInfinity(m_sigma)) throw new ArgumentException(string.Format("Sigma argument must be a finite number > 0 (got {0:G}).", m_sigma));
        }

        public override bool discrete() { return false; }

        public override bool tail_left() { return false; }

        public double sigma() { return m_sigma; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return range();
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            double sigmasqr = m_sigma * m_sigma;
            return x * (Math.Exp(-(x * x) / (2 * sigmasqr))) / sigmasqr;
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : 0;
            if (RHS) return find_pdf_inv(p, mode(), quantilec(double.Epsilon), false);
            return find_pdf_inv(p, 0, mode(), true);
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            return -XMath.expm1(-x * x / (2 * m_sigma * m_sigma));
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            return Math.Exp(-x * x / (2 * m_sigma * m_sigma));
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 0) return 0;
            if (p == 1) throw new OverflowException();
            return Math.Sqrt(-2 * m_sigma * m_sigma * XMath.log1p(-p));
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            if (q == 1) return 0;
            if (q == 0) throw new OverflowException();
            return Math.Sqrt(-2 * m_sigma * m_sigma * Math.Log(q));
        }

        public override double mean()
        {
            return m_sigma * Math.Sqrt(Math.PI / 2);
        }

        public override double variance()
        {
            return (4 - Math.PI) * m_sigma * m_sigma / 2;
        }

        public override double mode()
        {
            return m_sigma;
        }

        public override double median()
        {
            return Math.Sqrt(Math.Log(4.0)) * m_sigma;
        }

        public override double skewness()
        {
            return 0.63111065781893713819189935154422777984404221106391;
        }

        public override double kurtosis()
        {
            return 3.2450893006876380628486604106197544154170667057995;
        }

        public override double kurtosis_excess()
        {
            return 0.2450893006876380628486604106197544154170667057995;
        }
    }
}
