using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class extreme_value_distribution : distribution
    {
        double m_a, m_b;

        public extreme_value_distribution(double a, double b)
        {
            m_a = a;
            m_b = b;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (double.IsInfinity(m_a)) throw new ArgumentException(string.Format("Location argument must be a finite number (got {0:G}).", m_a));
            if (m_b <= 0 || double.IsInfinity(m_b)) throw new ArgumentException(string.Format("Scale argument must be a finite number > 0 (got {0:G}).", m_b));
        }


        public override bool discrete() { return false; }

        public double location() { return m_a; }

        public double scale() { return m_b; }

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
            return Math.Exp((m_a - x) / m_b) * Math.Exp(-Math.Exp((m_a - x) / m_b)) / m_b;
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : -double.MaxValue;
            if (RHS) return find_pdf_inv(p, mode(), quantilec(double.Epsilon), false);
            else return find_pdf_inv(p, quantile(double.Epsilon), mode(), true);
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            return Math.Exp(-Math.Exp((m_a - x) / m_b));
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            return -XMath.expm1(-Math.Exp((m_a - x) / m_b));
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 0 || p == 1) throw new OverflowException();
            return m_a - Math.Log(-Math.Log(p)) * m_b;
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            if (q == 0 || q == 1) throw new OverflowException();
            return m_a - Math.Log(-XMath.log1p(-q)) * m_b;
        }

        public override double mean()
        {
            return m_a + XMath.euler * m_b;
        }

        //Variance supplied by base class
        public override double standard_deviation()
        {
            return Math.PI * m_b / Math.Sqrt(6.0);
        }

        public override double mode()
        {
            return m_a;
        }

        public override double median()
        {
            return m_a - m_b * XMath.ln_ln_two;
        }

        public override double skewness()
        {
            return 1.1395470994046486574927930193898461120875997958366;
        }

        public override double kurtosis()
        {
            return 27.0 / 5.0;
        }

        public override double kurtosis_excess()
        {
            return 12.0 / 5.0;
        }

    }
}
