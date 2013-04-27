using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class pareto_distribution : distribution
    {
        double m_scale, m_shape;

        public pareto_distribution(double shape, double scale)
        {
            m_scale = scale;
            m_shape = shape;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_shape <= 0 || double.IsInfinity(m_shape)) throw new ArgumentException(string.Format("Shape argument must be a finite number > 0 (got {0:G}).", m_shape));
            if (m_scale <= 0 || double.IsInfinity(m_scale)) throw new ArgumentException(string.Format("Scale argument must be a finite number > 0 (got {0:G}).", m_scale));
        }

        public double scale() { return m_scale; }

        public double shape() { return m_shape; }

        public override bool discrete() { return false; }

        public override bool LHS() { return false; }

        public override bool tail_left() { return false; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(m_scale, double.MaxValue);
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            if (x < m_scale) return 0;
            return m_shape * Math.Pow(m_scale, m_shape) / Math.Pow(x, m_shape + 1);
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            if (!RHS) throw new Exception("The Pareto distribution is one-sided (RHS only) - no LHS inverse PDF is available.");
            base.pdf_inv(p, RHS);
            if (p == 0) return double.MaxValue;
            double x = Math.Pow((m_shape * Math.Pow(m_scale, m_shape) / p), 1.0 / (m_shape + 1));
            if (double.IsInfinity(x)) x = support().v2;
            return x;
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if (x < m_scale) return 0;
            return -XMath.powm1(m_scale / x, m_shape);            
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (x < m_scale) return 1;
            return Math.Pow(m_scale / x, m_shape);
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 0) return m_scale;
            if (p == 1) return double.PositiveInfinity;
            return m_scale / (Math.Pow(1 - p, 1 / m_shape));
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            if (q == 1) return m_scale;
            if (q == 0) return double.PositiveInfinity;
            return m_scale / (Math.Pow(q, 1 / m_shape));
        }

        public override double mean()
        {
            if (m_shape > 1) return m_shape * m_scale / (m_shape - 1);
            return double.PositiveInfinity;
        }

        public override double mode()
        {
            return m_scale;
        }

        public override double median()
        {
            return m_scale * Math.Pow(2.0, 1 / m_shape);
        }

        public override double variance()
        {
            if (m_shape <= 2) throw new Exception(string.Format("Pareto distribution: Variance is not defined for m_shape <= 2 (m_shape = {0:G}).", m_shape));
            return (m_scale * m_scale * m_shape) / ((m_shape - 1) * (m_shape - 1) * (m_shape - 2));
        }

        public override double skewness()
        {
            if (m_shape <= 3) throw new Exception(string.Format("Pareto distribution: Skewness is not defined for m_shape <= 3 (m_shape = {0:G}).", m_shape));
            return Math.Sqrt((m_shape - 2) / m_shape) * 2 * (m_shape + 1) / (m_shape - 3);
        }

        public override double kurtosis()
        {
            if (m_shape <= 4) throw new Exception(string.Format("Pareto distribution: Kurtosis is not defined for m_shape <= 4 (m_shape = {0:G}).", m_shape));
            return 3 * ((m_shape - 2) * (3 * m_shape * m_shape + m_shape + 2)) / (m_shape * (m_shape - 3) * (m_shape - 4));
        }

        public override double kurtosis_excess()
        {
            if (m_shape <= 4) throw new Exception(string.Format("Pareto distribution: Kurtosis is not defined for m_shape <= 4 (m_shape = {0:G}).", m_shape));
            return 6 * ((m_shape * m_shape * m_shape) + (m_shape * m_shape) - 6 * m_shape - 2) / (m_shape * (m_shape - 3) * (m_shape - 4));
        }
    
    }
}
