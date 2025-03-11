using System;
using System.Numerics;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.IO; // ��� ������ � �������

class Program : Form
{
    private Chart chart;

    static void Main()
    {
        Application.EnableVisualStyles();
        Application.Run(new Program());
    }

    public Program()
    {
        // ������� ����� � ������
        this.Text = "������� �������";
        this.Size = new System.Drawing.Size(800, 600);

        chart = new Chart();
        chart.Dock = DockStyle.Fill;

        var chartArea = new ChartArea();
        chart.ChartAreas.Add(chartArea);

        var series1 = new Series("f(x) = |sin(x)|");
        var series2 = new Series("p*(x) (���)");
        var series3 = new Series("p*(x) (������)");

        series1.ChartType = SeriesChartType.Line;
        series2.ChartType = SeriesChartType.Line;
        series3.ChartType = SeriesChartType.Line;

        series1.Color = Color.Red;
        series2.Color = Color.Blue;
        series3.Color = Color.Magenta;

        chart.Series.Add(series1);
        chart.Series.Add(series2);
        chart.Series.Add(series3);

        this.Controls.Add(chart);

        // ���������
        int n = 16; // ������� ������������������� ����������
        int N = n * 2; // ���������� ����� ��� ������������� 
        double step = 2 * Math.PI / N; // ��� �������������

        // ������������� ������� f(x) = |sin(x)|
        Complex[] samples = new Complex[N];
        for (int i = 0; i < N; i++)
        {
            double x = i * step;
            samples[i] = Math.Abs(Math.Sin(x));
        }

        // ���������� ���
        Complex[] fftCoefficients = FFT(samples);

        // ���������� ������ n ������������� �����
        double[] realCoefficients = new double[n];
        for (int k = 0; k < n; k++)
        {
            realCoefficients[k] = fftCoefficients[k].Real / N;
        }

        // ��������� ���� ��� ������
        using (StreamWriter writer = new StreamWriter("output.txt"))
        {
            // ���������� ��������� � ����
            writer.WriteLine("x\tf(x) = |sin(x)|\tp*(x) (���)\tp*(x) (������)");

            // ���������� �������� f(x), p*(x) � ������� ������� �� ����������� �����
            for (int i = 0; i <= N; i++)
            {
                double x = i * step;
                double fx = Math.Abs(Math.Sin(x)); // �������� f(x)
                double px_fft = ApproximationPolynomial(realCoefficients, x); // �������� p*(x) (���)
                double px_exact = ExactSolution(x); // �������� p*(x) (������)

                // ��������� ����� � �������
                chart.Series["f(x) = |sin(x)|"].Points.AddXY(x, fx);
                chart.Series["p*(x) (���)"].Points.AddXY(x, px_fft);
                chart.Series["p*(x) (������)"].Points.AddXY(x, px_exact);

                // ���������� �������� � ����
                writer.WriteLine($"{x:F4}\t{fx:F6}\t{px_fft:F6}\t{px_exact:F6}");
            }

            // ������������� ���������� ������������ ���
            writer.WriteLine("\n������������ ���:");
            writer.WriteLine("k\tRe(FFT[k])\tIm(FFT[k])");
            for (int k = 0; k < fftCoefficients.Length; k++)
            {
                writer.WriteLine($"{k}\t{fftCoefficients[k].Real:F6}\t{fftCoefficients[k].Imaginary:F6}");
            }
        }

        Console.WriteLine("������ ������� �������� � ���� output.txt");
    }

    // ���������� �������� �������������� ����� (���)
    static Complex[] FFT(Complex[] x)
    {
        int N = x.Length;
        if (N <= 1)
            return x;

        // ���������� �� ������ � �������� ��������
        Complex[] even = new Complex[N / 2];
        Complex[] odd = new Complex[N / 2];
        for (int i = 0; i < N / 2; i++)
        {
            even[i] = x[2 * i];
            odd[i] = x[2 * i + 1];
        }

        // ����������� ����� ���
        Complex[] fftEven = FFT(even);
        Complex[] fftOdd = FFT(odd);

        // ����������� �����������
        Complex[] fftResult = new Complex[N];
        for (int k = 0; k < N / 2; k++)
        {
            double angle = -2 * Math.PI * k / N;
            Complex t = Complex.FromPolarCoordinates(1, angle) * fftOdd[k];
            fftResult[k] = fftEven[k] + t;
            fftResult[k + N / 2] = fftEven[k] - t;
        }

        return fftResult;
    }

    // ������� ��� ���������� �������� ���������� ���������� ����������� p*(x) (���)
    static double ApproximationPolynomial(double[] coefficients, double x)
    {
        double result = coefficients[0]; // ������� �����������
        for (int k = 1; k < coefficients.Length; k++)
        {
            result += coefficients[k] * Math.Cos(k * x); // ��������� ���������
        }
        return result;
    }

    // ������� ��� ���������� ������� ������� p*(x)
    static double ExactSolution(double x)
    {
        double result = 2 / Math.PI; // ����������� ����
        for (int k = 1; k <= 7; k++)
        {
            result -= (4 / Math.PI) * Math.Cos(2 * k * x) / (4 * k * k - 1);
        }
        return result;
    }
}