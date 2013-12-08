import jbios.analysis.AnalysisResult;
import jbios.analysis.AnalysisResult.ResultSeriesGroup;
import jbios.plugin.AnalysisMethod;
import jbios.plugin.InputParametersDialog;
import jbios.plugin.InputParametersException;
import jbios.signal.Signal;

public class FFT_ implements AnalysisMethod {
@Override
    public String getDescription() {
        return "Fast Fourier Transform";
    }
@Override    
    public void inputParameters() throws InputParametersException {
    }
@Override
    public AnalysisResult doAnalysis(Signal[] signals) throws InterruptedException {
        
        double[] amostras = signals[0].getSamples(); /* extrai os valores do sinal e atribiu a um vetor double */
        
        double[] signal = PotenciaDeDois(amostras); /* verifica se o vetor possui 2^x valores, corrige caso nao tenha */
        
        double[] imaginary = new double[signal.length]; /* cria um novo vetor double para a frequência */
        
        for(int i = 0; i < signal.length; i++){
            imaginary[i] = 0; /* preenche o vetor */
        }
        
        // Aplica a transformada rápida de Fourier
        double[][] FFTresult = new double [2][signal.length];
        
        FFTresult = FFT(signal, imaginary); /* aplica a tranformada rápida de fourier */
        
        double[] frequencia = new double[signal.length]; /* cria um novo vetor double para a frequência */
        
        double Hz = signals[0].getSamplingRate();
        
        for(int i = 0; i<signal.length; i++){
            frequencia[i] = i * (Hz / signal.length); /* preenche o vetor */
        }
        
        // Calcula os valores de magnitude e fase da transformada de Fourier
        double[] magnitude = CalculoMagnitude(FFTresult);
        double[] fase = CalculoFase(FFTresult);
        
        // Calcula a potência para cada frequência de banda
        double[] potencias = Potencias_Banda(signal, frequencia, Hz);
        
        AnalysisResult resultFFT = new AnalysisResult("Plugin FFT"); /* define o nome da aba de resulatados */
        
        // Criando grupo de series.
        // Grupo de resultados da parte real
        ResultSeriesGroup group = resultFFT.new ResultSeriesGroup(ResultSeriesGroup.XY,
                                                                   "Parte Real", 
                                                                   "Frequência (Hz)", 
                                                                   "Valor", 
                                                                   0);
        
        // Grupo de resultados da parte imaginária
        ResultSeriesGroup group2 = resultFFT.new ResultSeriesGroup(ResultSeriesGroup.XY,
                                                                   "Parte Imaginária", 
                                                                   "Frequência (Hz)", 
                                                                   "Valor", 
                                                                   0);
        
        // Grupo de resultados de magnitude x fase
        ResultSeriesGroup group3 = resultFFT.new ResultSeriesGroup(ResultSeriesGroup.XY,
                                                                   "Magnitude x Fase", 
                                                                   "Fase", 
                                                                   "Magnitude", 
                                                                   0);
        
        // Adicionando uma série ao grupo.
        group.addResultSeries("Original Signal", frequencia, FFTresult[0], null, null, false, null);
        
        group2.addResultSeries("Original Signal", frequencia, FFTresult[1], null, null, false, null);
        
        group3.addResultSeries("Original Signal", fase, magnitude, null, null, false, null);
        
        // Adicionando o grupo aos resultados.
        resultFFT.addResultSeriesGroup(group);
        resultFFT.addResultSeriesGroup(group2);
        resultFFT.addResultSeriesGroup(group3);
        
        // Adicionando resultados de potência média
        resultFFT.addResultSingleValue("Potência média para banda de frequência Beta", "-", potencias[0]);
        resultFFT.addResultSingleValue("Potência média para banda de frequência Alfa", "-", potencias[1]);
        resultFFT.addResultSingleValue("Potência média para banda de frequência Theta", "-", potencias[2]);
        resultFFT.addResultSingleValue("Potência média para banda de frequência Delta", "-", potencias[3]);
        
        // É recomendado o uso de Thread.sleep(0) antes de retornar o resultado.
        Thread.sleep(0);
        return resultFFT;
	
    }
    
    // Aplica a transformada rápida de Fourier de Cooley e Tukey.
    public static double[][] FFT(double real[], double imaginary[]) {
        int shift = 0;
        int size = real.length;
        int m = (int)Math.floor(Math.log(size) / Math.log(2D));
        int n = 1 << m;
        double Imarg[] = new double[n];
        double Rearg[] = new double[n];
        double arg0 = 2 * Math.PI / (double)n;
        double[][] FFTresult = new double [2][real.length];
        for(int i = 0; i < n; i++)
        {
            double arg = arg0 * (double)i;
            Rearg[i] = Math.cos(arg);
            Imarg[i] = -Math.sin(arg);
        }

        int j;
        for(int i = j = shift; i < (shift + n) - 1; i++)
        {
            if(i < j)
            {
                double Retmp = real[i];
                double Imtmp = imaginary[i];
                real[i] = real[j];
                imaginary[i] = imaginary[j];
                real[j] = Retmp;
                imaginary[j] = Imtmp;
            }
            int k;
            for(k = n >> 1; k + shift <= j; k /= 2)
                j -= k;

            j += k;
        }

        int stepsize = 1;
        for(int shifter = m - 1; stepsize < n; shifter--) {
            for(j = shift; j < shift + n; j += stepsize << 1) {
                for(int i = 0; i < stepsize; i++) {
                    int i_j = i + j;
                    int i_j_s = i_j + stepsize;
                    double Retmp;
                    if(i > 0)
                    {
                        Retmp = Rearg[i << shifter] * real[i_j_s] - Imarg[i << shifter] * imaginary[i_j_s];
                        imaginary[i_j_s] = Rearg[i << shifter] * imaginary[i_j_s] + Imarg[i << shifter] * real[i_j_s];
                        real[i_j_s] = Retmp;
                    }
                    Retmp = real[i_j] - real[i_j_s];
                    double Imtmp = imaginary[i_j] - imaginary[i_j_s];
                    real[i_j] += real[i_j_s];
                    imaginary[i_j] += imaginary[i_j_s];
                    real[i_j_s] = Retmp;
                    imaginary[i_j_s] = Imtmp;
                }

            }

            stepsize <<= 1;
        }
        FFTresult[0] = real;
        FFTresult[1] = imaginary;
        return(FFTresult);
    }
    
    // Verifica se o sinal possui N = 2^x pontos, completa com zeros até o próximo N que satisfaça a condição.
    public static double[] PotenciaDeDois(double signal[]){
        int x = 1;
        if(Integer.bitCount(signal.length) == 1){ /* verifica se o sinal possui n = 2^x */
            return(signal);
        }
        else{
            while(x < signal.length){
                    x = x*2;
            }
            double[] pow_2 = new double[x]; /* cria um vetor com n=x^2 e copia os dados para o mesmo */
            System.arraycopy(signal, 0, pow_2, 0, signal.length);
            for(int j=signal.length; j < x; j++){
                pow_2[j] = 0; /* completa o vetor com zeros */
            }
        return(pow_2);
        }
    }
    
    // Calcula o espectro de magnitude da transformada.
    public static double[] CalculoMagnitude(double FFTsignal[][]) {
        double[] magnitude = new double[FFTsignal[0].length];
        for(int i=0; i < FFTsignal[0].length; i++){
            magnitude[i] = Math.sqrt( Math.pow(FFTsignal[0][i], 2) + Math.pow(FFTsignal[1][i], 2) );
        }
        return(magnitude);
    }
    
    // Calcula o espectro de fase da transformada.
    public static double[] CalculoFase(double FFTsignal[][]) {
        double[] fase = new double[FFTsignal[0].length];
        for(int i=0; i < FFTsignal[0].length; i++){
            fase[i] = Math.atan( FFTsignal[1][i] / FFTsignal[0][i] );
        }
        return(fase);
    }

    // Calcula os pontos de inicio e fim das bandas de frequencia.
    public static double[] Potencias_Banda(double signal[], double frequencia[], double Hz) {
        int[] inicio_fim = new int[8];
        double[] resultados = new double[4];
        // ondas Beta
        inicio_fim[0] = (int)Math.floor( (14 * frequencia.length) / Hz );
        inicio_fim[1] = (int)Math.floor( (30 * frequencia.length) / Hz );
        // ondas Alfa
        inicio_fim[2] = (int)Math.floor( (8 * frequencia.length) / Hz );
        inicio_fim[3] = (int)Math.floor( (13 * frequencia.length) / Hz );
        // ondas Theta
        inicio_fim[4] = (int)Math.floor( (4 * frequencia.length) / Hz );
        inicio_fim[5] = (int)Math.floor( (7 * frequencia.length) / Hz );
        // ondas Delta
        inicio_fim[6] = (int)Math.floor( (0.5 * frequencia.length) / Hz );
        inicio_fim[7] = (int)Math.floor( (3 * frequencia.length) / Hz );

        for(int j=0; j <= 6; j=j+2){
            resultados[j/2] = Potencia_Media(inicio_fim[j], inicio_fim[j+1], signal);
        }
        return(resultados);
    }
    
    //Calcula a potencia média de um sinal.
    public static double Potencia_Media(int inicio, int fim, double signal[]){
        double potencia = 0;
        for(int i=inicio; i <= fim; i++) {
            potencia = potencia + ( Math.pow(signal[i], 2) / ( 2 * (fim - inicio) ) );
        }
        return(potencia);
    }
}