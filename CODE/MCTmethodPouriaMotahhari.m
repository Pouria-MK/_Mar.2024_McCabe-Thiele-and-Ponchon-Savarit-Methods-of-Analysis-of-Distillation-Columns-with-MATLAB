clc
clear
close all
syms x
%% Pouria Motahhari 
%Unit Operations Project

%#######################################################################
%#################Data Inputs and Error Handlings Part #################

%feed1:
F1=input('F1= ');
z1=input('z1= ');
q1=input('q1= ');

%feed2:
z2=z1;
q2=q1;
F2=input('Double-Feed (OPTIONAL): F2= ');
if F2>0
    z2=input('z2= ');
    if z2==z1
        error("z1 equals z2... then it's a single feed distillation!")
    end

    q2=input('q2= ');
end

if F2<0 || F1<0 || F1==0
    error("F1 must be must be higher than 0 and F2 must be must be higher or equal to 0!")
end
if z1<=0 || z1>=1 || z2<0 || z1>=1
    error("z1 & z2 must be between 0 to 1 and z1 must be positive!")
end

%in distillation columns: z1>z2
if z2>z1
    swap(z1,z2)
    swap(q1,q2)
    swap(F1,F2)
end

%products
xD=input('xD= ');
xB=input('xB= ');

if xD<=0 || xD>=1 || xB<=0 || xB>=1
    error("xD & xB Domains is (0,1)")
end

%compositions order:
if ((xB<z2) && (z2<=z1) && (z1<xD))==0
    error("Compositions order must be: xB<zF2<=zF1<xD")
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

%equilibrium curve fitting:
[eq_curve,gof_a]=fit(xe,ye,'poly5'); %gof= goodness of fit

eq_curve_sym=subs(str2sym(formula(eq_curve)), ...
    coeffnames(eq_curve), ...
    num2cell(coeffvalues(eq_curve).')); %for fitted polynomial curve equation


%r_square value to measure the wellness of our fit model (MORE INFO ON REPORT)
%for equilibrium curve:
gof_poly_array=table2array(struct2table(gof_a));
r_square_poly=gof_poly_array(2); %if r^2>0.95: good fitting model for equilibrium curve





%#########################################################################
%################# Distillate and Bottoms Flow Rate Part #################
%Solve for D & B: (A_solve*X_solve=C_solve) "2 equations, 2 unknowns"
A_solve=[1 1
        xD xB];
C_solve=[F1+F2
        F1*z1+F2*z2];
X_solve=A_solve^-1*C_solve;

D=X_solve(1);
B=X_solve(2);
%D & B to be printed later




%######################################################################
%################# Feed-Lines(q-lines) Equations Part #################
feedline1=x*q1/(q1-1)-z1/(q1-1);
feedline2=x*q2/(q2-1)-z2/(q2-1);



%#########################################################
%################# Nmin Calculation Part #################
%Nmin calculation by assuming the rectifying section line is line x=y (Reflux-->inf)
min_stage_comp_x=[vpasolve(eq_curve_sym==xD,[0,1])]; %vpasolve(eq_curve_sym=xD,[0,1]) is the x value of first stage
min_stage_comp_y=[xD]; %y value of first stage = xD (or yD)

while min_stage_comp_x(end) > xB
    min_stage_comp_y(end+1)= min_stage_comp_x(end); %y value of next stage added to array y's
    min_stage_comp_x(end+1)= vpasolve( eq_curve_sym==min_stage_comp_y(end) ,[0,1]); %x value of next stage added to array x'es
end

Nmin=length( min_stage_comp_y ); %  (with reboiler)

%to calculate Nmin float part we can use this equation:
%lastx=(min_stage_comp_x(end-1)-xB)/(min_stage_comp_x(end-1)-min_stage_comp_x(end));





%########################################################
%#################Rmin Calculation Part #################
%Rmin calculation by feed-line1, equilibrium-curve and rectifying-section-line intersection:
if q1==1
    x_intersect_min= z1; %if q1 equals 1, feedline1 will be a vertical line
else
    x_intersect_min= vpasolve(feedline1==eq_curve_sym,[0,1]); %q-line and eqillibrium_curve x intersection
end
y_intersect_min= subs(eq_curve_sym,x_intersect_min);


%minimum rectifying section line slope= Rmin/(Rmin+1)= (yD-y_intersect)/(xD-x_intersect):
Rmin=double( (xD-y_intersect_min) / (y_intersect_min-x_intersect_min) ); %double command= symbolic to numeric


%====Correction for Rmin (if the distillation column runs on double feed)====
%========================== (MORE INFO ON REPORT) ==========================
if F2~=0

    %middle section operating line equation:
    %equation of this operating line is: y= x*(L'/G') + (D*xD-F1*zF1)/(G')     ;slope=L'/G'
    Lmin=Rmin*D;
    Gmin=Lmin+D;
    L_prime_min=q1*F1+Lmin;   % q1= (L'-L)/F1
    G_prime_min=Gmin+F1*(q1-1);  % 1-q1= (G-G')/F1
    
    mid_feeds_operating_line_min= x*(L_prime_min/G_prime_min) + (D*xD-F1*z1)/G_prime_min ;
    
    %minimum middle section line and feedline2 intersection coordinates:
    if q2==1
        x_intersect_min2= z2; %if q1 equals 1, feedline1 will be a vertical line
    else
        x_intersect_min2= double( vpasolve( feedline2==mid_feeds_operating_line_min ,[0,1]) ); %q-line and eqillibrium_curve x intersection
    end
    y_intersect_min2= double( subs(mid_feeds_operating_line_min,x_intersect_min2) );
    
    
    if y_intersect_min2> double( subs(eq_curve_sym,x_intersect_min2) )
    
        if q2==1
        x_intersect_min= z2; %if q1 equals 1, feedline1 will be a vertical line
        else
        x_intersect_min= vpasolve(feedline2==eq_curve_sym,[0,1]); %q-line and eqillibrium_curve x intersection
        end
        y_intersect_min= subs(eq_curve_sym,x_intersect_min);
    
        %now that we have the coordinates of feedline2 and equilibruim curve
        %intersect (rectifying section line), we can calculate new G_bar, L_bar , G_prime, L_prime (from q2),
        %G and L (from q1) in minimum reflux condition:
    
        %stripping section interpole = (B*xB/G_bar) (based on 2 points of (xB,xB) and
        %feedline2 and eq_curve intersect: (more info on report)
        G_bar_min= -(B*xB)/( xB- xB* (y_intersect_min-xB)/(x_intersect_min-xB) );
    
        %stripping section slope= L_bar/G_bar :
        L_bar_min= G_bar_min* (y_intersect_min-xB)/(x_intersect_min-xB); %unused variable (educational purposes)
    
        %new G_prime_min & L_prime_min calculations by q2 definition:
        G_prime_min= G_bar_min+ F2*(1-q2);
        L_prime_min= L_bar_min- F2*q2; %unused variable (educational purposes)
    
        %new G_min & L_min calculations by q1 definition:
        G_min= G_prime_min+ F1*(1-q1);
        L_min= L_prime_min- F1*q1; %unused variable (educational purposes)
    
        %new minimum reflux calculation by new G_min & L_min:
        Rmin= G_min/D -1;
    
    end
end

%================== End of correction for double-feed Rmin ==================




%########################################################
%################# Operating Lines Part #################
%-------------------rectifying section line equation---------------------
Reflux=Rmin*reflux_ratio;
rectifying=x*Reflux/(Reflux+1)+xD/(Reflux+1);

%rectifying section line and feedline1 intersection coordinates:
if q1==1
    x_intersect1= z1; %if q1 equals 1, feedline1 will be a vertical line
else
    x_intersect1= double( vpasolve(feedline1==rectifying,[0,1]) ); %q-line and eqillibrium_curve x intersection
end
y_intersect1= double( subs(rectifying,x_intersect1) );


%--------------------middle section operating line equation--------------------
%if its a single feed distillation, because q-line1 and q-line2 are the
%same and matched, middle sectione line is actually a point and doesn't get
%shown in diagram (MORE INFO ON REPORT)

%equation of this operating line is: y= x*(L'/G') + (D*xD-F1*zF1)/(G')     ;slope=L'/G'
L=Reflux*D;
G=L+D;
L_prime=q1*F1+L;   % q1= (L'-L)/F1
G_prime=G+F1*(q1-1);  % 1-q1= (G-G')/F1
    
mid_feeds_operating_line= x*(L_prime/G_prime) + (D*xD-F1*z1)/G_prime ;
    
%middle section line and feedline2 intersection coordinates:
if q2==1
        x_intersect2= z2; %if q2 equals 1, feedline2 will be a vertical line
else
        x_intersect2= double( vpasolve( feedline2==mid_feeds_operating_line ,[0,1]) ); %q-line and eqillibrium_curve x intersection
end
y_intersect2= double( subs(mid_feeds_operating_line,x_intersect2) );



%----------------------stripping section line equation--------------------------
%from the coordination of 2 points (xB,xB) & (x_intersect2,y_intersect2)
stripping = (x-xB)*(y_intersect2-xB)/(x_intersect2-xB) +xB;



%########################################################
%################# NTP Calculation Part #################
%theoretical stages calculation by adding stage compositions in a while loop

%rectifying section:
stage_comp_x=[ vpasolve(eq_curve_sym==xD,[0,1]) ]; %vpasolve(eq_curve_sym=xD,[0,1]) is the x value of first stage
stage_comp_y=[xD]; %y value of first stage = xD (or yD)

%stages drawer array to plot: (FULL EXPLAINATION ON REPORT)
stage_drawer_x=[xD, stage_comp_x];
stage_drawer_y=[xD, stage_comp_y];


while stage_comp_x(end) > x_intersect1
    stage_comp_y(end+1)= double( subs( rectifying,stage_comp_x(end) ) ); %y value of next stage added to array y's
    stage_comp_x(end+1)= vpasolve( eq_curve_sym==stage_comp_y(end) ,[0,1]); %x value of next stage added to array x'es

    stage_drawer_x=[stage_drawer_x, stage_comp_x(end-1), stage_comp_x(end)];
    stage_drawer_y=[stage_drawer_y, stage_comp_y(end), stage_comp_y(end)];
end

feed_stage1= length( stage_comp_y ); %feed stage for feed 1


%middle (between 2 feeds) section: (if it's a single feed distillation, x_intersect1=x_intersect2 thus this part skips)
while stage_comp_x(end) > x_intersect2
    stage_comp_y(end+1)= double( subs( mid_feeds_operating_line,stage_comp_x(end) ) ); %y value of next stage added to array y's
    stage_comp_x(end+1)= vpasolve( eq_curve_sym==stage_comp_y(end) ,[0,1]); %x value of next stage added to array x'es

    stage_drawer_x=[stage_drawer_x, stage_comp_x(end-1), stage_comp_x(end)];
    stage_drawer_y=[stage_drawer_y, stage_comp_y(end), stage_comp_y(end)];
end

feed_stage2= length( stage_comp_y ); %feed stage for feed 2


%stripping section:
while stage_comp_x(end) > xB
    stage_comp_y(end+1)= double( subs( stripping,stage_comp_x(end) ) ); %y value of next stage added to array y's
    stage_comp_x(end+1)= vpasolve( eq_curve_sym==stage_comp_y(end) ,[0,1]); %x value of next stage added to array x'es

    stage_drawer_x=[stage_drawer_x, stage_comp_x(end-1), stage_comp_x(end)];
    stage_drawer_y=[stage_drawer_y, stage_comp_y(end), stage_comp_y(end)];
end


%draw last vertical line for reboiler:
stage_drawer_x=[stage_drawer_x, stage_drawer_x(end)];
stage_drawer_y=[stage_drawer_y, stage_drawer_x(end)];

%Number of theoretical stages:
NTP=length( stage_comp_y );




%################################################
%################# OUTPUTS PART #################
%diagram:
hold on
title("McCabe-Thiele Method Output Diagram")
xlabel('x= Mole fraction of most volatile component in liquid phase') 
ylabel('y= Mole fraction of most volatile component in vapor phase') 
box on
pbaspect([1,1,1]) %to make square diagram

fplot(eq_curve_sym,[0,1], "B") %equilibrium curve

plot(stage_drawer_x,stage_drawer_y,"r") %stages
text(stage_comp_x(1:end-1) -0.01 ,stage_comp_y(1:end-1) +0.01, string(1:numel(stage_comp_x)-1), 'FontSize',8 )  %stages numbering 
text(stage_comp_x(end) ,stage_comp_y(end)-0.005 , "reboiler" ,'FontSize',6) %reboiler marker

plot([z1 x_intersect1] ,[z1 y_intersect1],"k") %feed-line1
plot([z2 x_intersect2] ,[z2 y_intersect2],"k") %feed-line2

fplot(rectifying,[x_intersect1,xD],"g") %rectifying section line
fplot(mid_feeds_operating_line,[x_intersect2,x_intersect1],"g") %middle section line
fplot(stripping,[xB,x_intersect2],"g") %stripping section line
plot(xe,xe,"y","linewidth",1) %x=y line

legend({'Equilibrium curve','Stages','Feed-line','','','Operating Lines','','x=y'},'Location','southeast')
axis([0 1 0 1])
hold off

%numerical answers and stats;
disp( "Distillate: D=  " +num2str(D) )
disp( "Bottoms: B=  " +num2str(B) )
disp("--------------------")
disp( "Minimum Stages (or infinite reflux): Nmin=  " +num2str(Nmin) )
disp( "Minimum Reflux: Rmin=  " +num2str( double(Rmin) ) )
disp( "Number of Theoretical Equilibrium stages (with Reboiler): NTP=  " +num2str(NTP) )

if z1==z2
    disp( "Feed stage: Nf=  " +num2str(feed_stage1) )
else 
    disp( "Feed stage 1: Nf1=  " +num2str(feed_stage1) )
    disp( "Feed stage 2: Nf2=  " +num2str(feed_stage2) )    
end

disp("--------------------")
disp( "stage compositions (last stage represents the reboiler) :" )
stage_comp_x=double(stage_comp_x); %sym to double to display
stage_comp_x=permute( stage_comp_x, [2,1] ); %to make the array vertical
stage_comp_y=permute( stage_comp_y, [2,1] ); %to make the array vertical
stage_counter=permute( 1:numel(stage_comp_x), [2,1] );

disp("Stage "+ num2str(stage_counter)+ ...
    ":  x="+ num2str(stage_comp_x)+ ...
    " | y="+ num2str(stage_comp_y)) %stage compositions
