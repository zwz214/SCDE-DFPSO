 load fbias_data
 y=f_bias;
 y1=(-1400:100:1400);
for i=1:25
       str1='Sheet';
       str2=num2str(Ni);
       str=strcat(str1,str2);
       de_trend=xlsread('fbest_scde_cec2013_2.xlsx',str,'A:A');
       fdr_trend=xlsread('fbest_scde_cec2013_2.xlsx',str,'D:D');
       de_trend=de_trend+y(i)-y1(i);
       fdr_trend=fdr_trend+y(i)-y1(i);
       xlswrite('fbest_scde_cec2013_2.xlsx',de_trend,str,'A1')
       xlswrite('fbest_scde_cec2013_2.xlsx',fdr_trend,str,'D1')
       de=xlsread('fbest_scde_cec2013_val_2.xlsx',str,'A:A');
       fdr=xlsread('fbest_scde_cec2013_val_2.xlsx',str,'D:D');
       de=de+y(i)-y1(i);
       fdr=fdr+y(i)-y1(i);
       xlswrite('fbest_scde_cec2013_val_2.xlsx',de,str,'a1')
       xlswrite('fbest_scde_cec2013_val_2.xlsx',fdr,str,'d1')
end