for i=1:10
       Ni=[1,3,5,6,11,12,18,22,25,27];
       str1='Sheet';
       str2=num2str(Ni(i));
       str=strcat(str1,str2);
       de_tre=xlsread('cec2013_trend_9',str,'B1£ºB3000');
       epso_tre=xlsread('cec2013_trend_9',str,'C1:C3000');
       fdr_tre=xlsread('cec2013_trend_9',str,'D1:D3000');
       DeDfpso_tre=xlsread('cec2013_trend_9',str,'A1:A3000');
       code_tre=xlsread('cec2013_trend_9',str,'E1:E3000');
       scdedfpso_tre=xlsread('cec2013_trend_9',str,'F1:F3000');
       x1=xlsread('cec2013_trend_9',str,'G1:G3000');
       x=floor(linspace(1,3000,50));
       x2=x1(x,:);
       de=de_tre(x,:);
       epso=epso_tre(x,:);
       fdr=fdr_tre(x,:);
       dedfpso=DeDfpso_tre(x,:);
       code=code_tre(x,:);
       scdedfpso=scdedfpso_tre(x,:);
       figure('position',[50,50,500,450])
       plot(x2,de,'b--^')
       hold on
       plot(x2,epso,'b--v')
       hold on
       plot(x2,fdr,'b--p')
       hold on
       plot(x2,dedfpso,'b--o')
       hold on
       plot(x2,code,'b--s')
       hold on
       plot(x2,scdedfpso,'b--d')
%        grid on
       xlabel('FES');
       ylabel('f(x)-f(x*)');
       legend('DE','EPSO','FDR','DE-DFPSO','CODE','SCDE-DFPSO')
end