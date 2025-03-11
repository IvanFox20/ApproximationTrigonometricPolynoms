using System;
using System.Numerics;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.IO; // Для работы с файлами

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
        // Создаем форму и график
        this.Text = "Графики функций";
        this.Size = new System.Drawing.Size(800, 600);

        chart = new Chart();
        chart.Dock = DockStyle.Fill;

        var chartArea = new ChartArea();
        chart.ChartAreas.Add(chartArea);

        var series1 = new Series("f(x) = |sin(x)|");
        var series2 = new Series("p*(x) (БПФ)");
        var series3 = new Series("p*(x) (точное)");

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

        // Параметры
        int n = 16; // Степень тригонометрического многочлена
        int N = n * 2; // Количество точек для дискретизации 
        double step = 2 * Math.PI / N; // Шаг дискретизации

        // Дискретизация функции f(x) = |sin(x)|
        Complex[] samples = new Complex[N];
        for (int i = 0; i < N; i++)
        {
            double x = i * step;
            samples[i] = Math.Abs(Math.Sin(x));
        }

        // Применение БПФ
        Complex[] fftCoefficients = FFT(samples);

        // Извлечение первых n коэффициентов Фурье
        double[] realCoefficients = new double[n];
        for (int k = 0; k < n; k++)
        {
            realCoefficients[k] = fftCoefficients[k].Real / N;
        }

        // Открываем файл для записи
        using (StreamWriter writer = new StreamWriter("output.txt"))
        {
            // Записываем заголовок в файл
            writer.WriteLine("x\tf(x) = |sin(x)|\tp*(x) (БПФ)\tp*(x) (точное)");

            // Вычисление значений f(x), p*(x) и точного решения на равномерной сетке
            for (int i = 0; i <= N; i++)
            {
                double x = i * step;
                double fx = Math.Abs(Math.Sin(x)); // Значение f(x)
                double px_fft = ApproximationPolynomial(realCoefficients, x); // Значение p*(x) (БПФ)
                double px_exact = ExactSolution(x); // Значение p*(x) (точное)

                // Добавляем точки в графики
                chart.Series["f(x) = |sin(x)|"].Points.AddXY(x, fx);
                chart.Series["p*(x) (БПФ)"].Points.AddXY(x, px_fft);
                chart.Series["p*(x) (точное)"].Points.AddXY(x, px_exact);

                // Записываем значения в файл
                writer.WriteLine($"{x:F4}\t{fx:F6}\t{px_fft:F6}\t{px_exact:F6}");
            }

            // Дополнительно записываем коэффициенты БПФ
            writer.WriteLine("\nКоэффициенты БПФ:");
            writer.WriteLine("k\tRe(FFT[k])\tIm(FFT[k])");
            for (int k = 0; k < fftCoefficients.Length; k++)
            {
                writer.WriteLine($"{k}\t{fftCoefficients[k].Real:F6}\t{fftCoefficients[k].Imaginary:F6}");
            }
        }

        Console.WriteLine("Данные успешно записаны в файл output.txt");
    }

    // Реализация быстрого преобразования Фурье (БПФ)
    static Complex[] FFT(Complex[] x)
    {
        int N = x.Length;
        if (N <= 1)
            return x;

        // Разделение на четные и нечетные элементы
        Complex[] even = new Complex[N / 2];
        Complex[] odd = new Complex[N / 2];
        for (int i = 0; i < N / 2; i++)
        {
            even[i] = x[2 * i];
            odd[i] = x[2 * i + 1];
        }

        // Рекурсивный вызов БПФ
        Complex[] fftEven = FFT(even);
        Complex[] fftOdd = FFT(odd);

        // Объединение результатов
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

    // Функция для вычисления значения многочлена наилучшего приближения p*(x) (БПФ)
    static double ApproximationPolynomial(double[] coefficients, double x)
    {
        double result = coefficients[0]; // Нулевой коэффициент
        for (int k = 1; k < coefficients.Length; k++)
        {
            result += coefficients[k] * Math.Cos(k * x); // Добавляем гармоники
        }
        return result;
    }

    // Функция для вычисления точного решения p*(x)
    static double ExactSolution(double x)
    {
        double result = 2 / Math.PI; // Константный член
        for (int k = 1; k <= 7; k++)
        {
            result -= (4 / Math.PI) * Math.Cos(2 * k * x) / (4 * k * k - 1);
        }
        return result;
    }
}