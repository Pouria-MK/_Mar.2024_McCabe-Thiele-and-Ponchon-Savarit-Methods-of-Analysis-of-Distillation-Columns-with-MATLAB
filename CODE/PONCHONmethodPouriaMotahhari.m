clc
clear
close all
syms x
%% Pouria Motahhari 

%Unit Operations Project

%#######################################################################
%#################Data Inputs and Error Handlings Part #################

%feed:
F=input('F= ');
z=input('z= ');
q=input('q= ');


%products:
xD=input('xD= ');
xB=input('xB= ');

if xD<=0 || xD>=1 || xB<=0 || xB>=1
    error("xD & xB Domains is (0,1)")
end

%compositions order:
if ((xB<z) && (z<xD))==0
    error("Compositions order must be: xB<zF<xD")
end

reflux_ratio=input('R/Rmin=  ');
if reflux_ratio<=1
    error("R/Rmin must be higher than 1 to operate!")
end




%###########################################################
%################# Equilibrium Curve Part #################
%Equilibrium Curve Data importation:
disp("Please upload the equilibrium file");
[file,path]=uigetfile('.txt');
%Path and File for my system= "C:\Users\p\Desktop\Unit Operations Project\Data & Properties\MeOH-Water.txt"
%Path and File string can be put in importdata command.

VLE=importdata( convertCharsToStrings([path file]) );
if isempty(VLE)==0
    disp("Equilibrium Data Loaded!")
    disp("=========================")
end

xe=VLE(:,1); %x_equilibrium
ye=VLE(:,2); %y_equilibrium
h_xe=VLE(:,3); %x_equilibrium
h_ye=VLE(:,4); %y_equilibrium

%equilibrium xy fitting:
[eq_xy,gof_a]=fit(xe,ye,'poly5'); %gof= goodness of fit

eq_xy_sym=subs(str2sym(formula(eq_xy)), ...
    coeffnames(eq_xy), ...
    num2cell(coeffvalues(eq_xy).')); %for fitted polynomial curve equation


%r_square value to measure the wellness of our fit model (MORE INFO ON REPORT)
%for equilibrium curve:
gof_poly_array=table2array(struct2table(gof_a));
r_square_poly=gof_poly_array(2); %if r^2>0.95: good fitting model for equilibrium curve



%equilibrium hx fitting:
[eq_hx,gof_b]=fit(xe,h_xe,'poly5'); %gof= goodness of fit

eq_hx_sym=subs(str2sym(formula(eq_hx)), ...
    coeffnames(eq_hx), ...
    num2cell(coeffvalues(eq_hx).')); %for fitted polynomial curve equation



%equilibrium hy fitting:
[eq_hy,gof_c]=fit(ye,h_ye,'poly5'); %gof= goodness of fit

eq_hy_sym=subs(str2sym(formula(eq_hy)), ...
    coeffnames(eq_hy), ...
    num2cell(coeffvalues(eq_hy).')); %for fitted polynomial curve equation





%#########################################################################
%################# Distillate and Bottoms Flow Rate Part #################
%Solve for D & B: (A_solve*X_solve=C_solve) "2 equations, 2 unknowns"
A_solve=[1 1
        xD xB];
C_solve=[F
        F*z];
X_solve=A_solve^-1*C_solve;

D=X_solve(1);
B=X_solve(2);
%D & B to be printed later




%#########################################################
%################# Nmin Calculation Part ################# 
% (same as mccabe-thiele method)
Nmin=1;
y_counter_min = xD;
x_counter_min = vpasolve( eq_xy_sym==y_counter_min ,[0,1]);

while x_counter_min > xB
    y_counter_min = x_counter_min;
    x_counter_min = vpasolve( eq_xy_sym==y_counter_min ,[0,1]);
    Nmin=Nmin+1;
end


%########################################################
%#################Rmin Calculation Part #################

x_intersect_min= z; %q equals 1, x_feed is on hx curve
y_intersect_min= subs(eq_xy_sym,x_intersect_min); %tie-line to x_intersect_min

hx_intersect_min= subs(eq_hx_sym,x_intersect_min); %enthalpy of x_intersect_min
hy_intersect_min= subs(eq_hy_sym,y_intersect_min); %enthalpy of y_intersect_min

feedline_min= ((hy_intersect_min-hx_intersect_min)/(y_intersect_min-x_intersect_min)) * (x-x_intersect_min) + hx_intersect_min;

%minimum reflux in ponchon-savarit method:
hy_d= subs(eq_hy_sym,xD);
hx_d= subs(eq_hx_sym,xD);
y_deltaD_min= subs(feedline_min,xD);

Rmin=double( (y_deltaD_min - hy_d) / (hy_d -hx_d) ); %double command= symbolic to numeric



%########################################################
%################### Feed Line Part #####################

Reflux = Rmin*reflux_ratio;
y_deltaD = (hy_d -hx_d)*Reflux + hy_d; %deltaD point

h_feed= subs(eq_hx_sym,z); %enthalpy of x_intersect_min (q=1)
feedline = ((y_deltaD-h_feed)/(xD-z)) * (x-z) + h_feed;

y_deltaW = subs(feedline,xB); %deltaW point


%########################################################
%########## NTP Calculation Part and plotting ###########
%diagram:
hold on
title("Ponchon-Savarit Method Output Diagram")
xlabel('x,y') 
ylabel('Enthalpy [J/Kmol]') 
box on

fplot(eq_hy_sym,[0,1], "B") %equilibrium curve enthalpy-y
fplot(eq_hx_sym,[0,1], "B") %equilibrium curve enthalpy-x
fplot(feedline,[xB,xD], "K") %feedline
xline(xB,"K") %xB
xline(xD,"K") %xB %xD


NTP=1;
y_counter = xD;
x_counter = vpasolve( eq_xy_sym==y_counter ,[0,1]);

%rectifying section
while x_counter>z
    
    hx_counter = subs(eq_hx_sym,x_counter);
    hy_counter= subs(eq_hy_sym,y_counter);
    plot([x_counter y_counter], [hx_counter hy_counter],"G")

    operating_line = ( ( y_deltaD-hx_counter ) /( xD-x_counter )) * (x-x_counter) + hx_counter;
    plot([x_counter xD],[hx_counter y_deltaD], "R") %operating-line

    y_counter = vpasolve( eq_hy_sym==operating_line ,[0,1]);
    x_counter = vpasolve( eq_xy_sym==y_counter ,[0,1]);

    NTP=NTP+1;
end


NF=NTP;


%stripping section
while x_counter>xB

    hx_counter = subs(eq_hx_sym,x_counter);
    hy_counter= subs(eq_hy_sym,y_counter);
    plot([x_counter y_counter], [hx_counter hy_counter],"G")

    operating_line = ( ( hx_counter-y_deltaW ) /( x_counter-xB )) * (x-xB) + y_deltaW;

    y_counter = vpasolve( eq_hy_sym==operating_line ,[0,1]);

    plot([xB y_counter],[y_deltaW hy_counter], "R") %operating-line

    x_counter = vpasolve( eq_xy_sym==y_counter ,[0,1]);

    NTP=NTP+1;
end

%reboiler stage (last stage)
hx_reb = subs(eq_hx_sym,x_counter);
hy_reb= subs(eq_hy_sym,y_counter);
plot([x_counter y_counter], [hx_reb hy_reb],"G")

operating_line = ( ( hx_reb-y_deltaW ) /( x_counter-xB )) * (x-xB) + y_deltaW;
plot([xB y_counter],[y_deltaW hy_reb], "R") %operating-line for reboiler

legend("Equillibrium curves",'',"Feedline",'','',"Tie-lines","Operating-lines",'','Location','southeast')
hold off

%numerical answers and stats;
disp( "Distillate: D=  " +num2str(D) )
disp( "Bottoms: B=  " +num2str(B) )
disp("--------------------")
disp( "Minimum Stages (or infinite reflux): Nmin=  " +num2str(Nmin) )
disp( "Minimum Reflux: Rmin=  " +num2str( double(Rmin) ) )
disp( "Number of Theoretical Equilibrium stages (with Reboiler): NTP=  " +num2str(NTP) )
disp( "Feed stage: Nf=  " +num2str(NF) )


