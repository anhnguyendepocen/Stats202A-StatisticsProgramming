*Final SAS project for STAT202A;
*Author: Anoosha Sagar;
*Complete the code;

*Import car data;
proc import out = work.data
datafile = "/folders/myfolders/regression_auto.csv"
dbms = csv replace; getnames=yes; datarow=2;
run;

*Compute the correlation between car length and mpg;
proc corr data=data;
var length weight;
run;

*Make a scatterplot of price (x-axis) and mpg (y-axis);
proc sgplot data=data;
title 'Scatter Plot for Price vs MPG';
scatter x = price y = mpg;
run;

*Make a box plot of mpg for foreign vs domestic cars;
proc sgplot data=data;
title 'MPG for Foreign and Domestic Cars';
vbox mpg / group=foreign;
run;

*Perform simple linear regression, y = mpg, x = price1;
* Do NOT include the intercept term;
proc reg data = data;
model mpg = price1 / noint;
run;

*Perform linear regression, y = mpg, x1 = length, x2 = length^2;
*Include the intercept term;
proc glm data=data;
model mpg = length length*length;
run;