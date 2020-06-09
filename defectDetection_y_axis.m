clc;close all;clear
%% 
index =1
str=dir(['C:\Users\lenovo\Desktop\HAM\Hum5\test_v_maskoutput',num2str(index,'%02d'),'b\*.jpg'])
[n,~]=size(str)
for i=1:n
    m1=imread(['C:\Users\lenovo\Desktop\HAM\Hum5\test_v_maskoutput',num2str(index,'%02d'),'b\p',num2str(i),'.jpg']); 
    BW1 = im2bw(m1);
    se = strel('disk',100);
    BW2 = imclose(BW1,se);
    s{i}=regionprops('table',BW2,{...
            'Centroid',...
            'MajorAxisLength',...
            'MinorAxisLength',...
            'Orientation'});
 end
 minor=[s{1} ; s{2}];
for j=3:n
   minor_n=[minor ; s{j}];
   minor=minor_n;
end
minor_new=table2array(minor);  
%% Add y_range to Table
u=1;
for i=1:n
    m1=imread(['C:\Users\lenovo\Desktop\HAM\Hum5\test_v_maskoutput',num2str(index,'%02d'),'b\p',num2str(i),'.jpg']); 
    BW1 = im2bw(m1);
    se = strel('disk',100);
    BW2 = imclose(BW1,se);
    s = regionprops(BW2,{...
            'Centroid',...
            'MajorAxisLength',...
            'MinorAxisLength',...
            'Orientation'});
    t = linspace(0,2*pi,50);
    % hold on
    for k = 1:length(s)
        for p = 1:45
        a = s(k).MajorAxisLength/2;
        b = s(k).MinorAxisLength/2;
        Xc = s(k).Centroid(1);
        Yc = s(k).Centroid(2);
        phi = deg2rad(-s(k).Orientation);
        x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
        y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);   
%         hold on
%         plot(x,y);
        end
        y1 = max(y);
        y2 = min(y);
        y_n{u}=y1-y2;
        y_max{u}=y1;
        y_min{u}=y2;
        u=u+1; 
    end
    y_nn=cell2mat(y_n);
    y_maxx=cell2mat(y_max);
    y_minn=cell2mat(y_min);
end 
y_range=reshape(y_nn,[u-1,1]);
y_range2=reshape(y_maxx,[u-1,1]);
y_range3=reshape(y_minn,[u-1,1]);

y_r=[minor_new,y_range];
y_r=[y_r,y_range2];
y_r=[y_r,y_range3];
%% delte extremely small ones of yrange
numRow=length(y_range);

m=max(y_range);
numSmall=numRow-n;
loop=0;
i=1;
while (loop~= numSmall)
    if y_range(i)<m/2
        y_r(i,:)=[];
        y_range(i)=[];
        y_range2(i)=[];
        y_range3(i)=[];
        loop=loop+1;
    else
        i=i+1;
    end
end 
%% Find out outliers of yrange
y_rangeNew=y_range;
Q1=quantile(y_rangeNew,0.25,1);
Q3=quantile(y_rangeNew,0.75,1);
IQR = iqr(y_rangeNew);
high=Q3+1.5*IQR;
low=Q3-1.5*IQR;
bound_up=Q3+3*IQR;
bound_low=Q1-3*IQR;  
outlierType1=y_rangeNew>high|y_rangeNew<low;
outlierType2=y_rangeNew>bound_up|y_rangeNew<bound_low;
y_rangeNew(outlierType1)=[];
while (any(outlierType2==1))
    Q1=quantile(y_rangeNew,0.25,1);
    Q3=quantile(y_rangeNew,0.75,1);
    IQR = iqr(y_rangeNew);
    top=Q3+1.5*IQR;
    low=Q3-1.5*IQR;
    bound_up=Q3+3*IQR;
    bound_low=Q1-3*IQR;
    outlierType1=y_rangeNew>top|y_rangeNew<low;
    outlierType2=y_rangeNew>bound_up|y_rangeNew<bound_low;
    y_rangeNew(outlierType1)=[];
end
%% find out outliers of ymin
y_rangeNew2=y_range2;
Q11=quantile(y_rangeNew2,0.25,1);
Q33=quantile(y_rangeNew2,0.75,1);
IQRR = iqr(y_rangeNew2);
highh=Q33+1.5*IQRR;
loww=Q33-1.5*IQRR;
bound_upp=Q33+3*IQRR;
bound_loww=Q11-3*IQRR;  
outlierType11=y_rangeNew2>highh|y_rangeNew2<loww;
outlierType22=y_rangeNew2>bound_upp|y_rangeNew2<bound_loww;
y_rangeNew2(outlierType11)=[];
while (any(outlierType22==1))
    Q11=quantile(y_rangeNew2,0.25,1);
    Q33=quantile(y_rangeNew2,0.75,1);
    IQRR = iqr(y_rangeNew2);
    topp=Q33+1.5*IQRR;
    loww=Q33-1.5*IQRR;
    bound_upp=Q33+3*IQRR;
    bound_loww=Q11-3*IQRR;
    outlierType11=y_rangeNew2>topp|y_rangeNew2<loww;
    outlierType22=y_rangeNew2>bound_upp|y_rangeNew2<bound_loww;
    y_rangeNew2(outlierType11)=[];
end
%% find out outliers of ymax
y_rangeNew3=y_range3;
Q111=quantile(y_rangeNew3,0.25,1);
Q333=quantile(y_rangeNew3,0.75,1);
IQRRR = iqr(y_rangeNew3);
highhh=Q333+1.5*IQRRR;
lowww=Q333-1.5*IQRRR;
bound_uppp=Q333+3*IQRRR;
bound_lowww=Q111-3*IQRRR;  
outlierType111=y_rangeNew3>highhh|y_rangeNew3<lowww;
outlierType222=y_rangeNew3>bound_uppp|y_rangeNew3<bound_lowww;
y_rangeNew3(outlierType111)=[];
while (any(outlierType222==1))
    Q111=quantile(y_rangeNew3,0.25,1);
    Q333=quantile(y_rangeNew3,0.75,1);
    IQRRR = iqr(y_rangeNew3);
    toppp=Q333+1.5*IQRRR;
    lowww=Q333-1.5*IQRRR;
    bound_uppp=Q333+3*IQRRR;
    bound_lowww=Q111-3*IQRRR;
    outlierType111=y_rangeNew3>toppp|y_rangeNew3<lowww;
    outlierType222=y_rangeNew3>bound_uppp|y_rangeNew3<bound_lowww;
    y_rangeNew3(outlierType111)=[];
end
%% save MinorAxisLength after delete all outliers
y_rangeT= array2table(y_rangeNew);
y_rangeT2= array2table(y_rangeNew2);
y_rangeT3= array2table(y_rangeNew3);
%% find the outliers indexs of yrange
j=1;
i=1;
indexx=0;
rowy_rangeNew=length(y_rangeNew);
diff=length(y_range)-length(y_rangeNew);
while (i~=length(y_range) || j<rowy_rangeNew)
    if (y_rangeNew(j)==y_range(i))
        j=j+1;
        i=i+1;
    else
        indexx=indexx+1;
        outlier(indexx)=i;
        i=i+1;
    end
end
if (indexx<diff)
    indexx=indexx+1;
    outlier(indexx)=length(y_range)
end
outlierT = array2table(outlier);
%% find the outliers indexs of ymin
jj=1;
ii=1;
indexxx=0;
rowy_rangeNew2=length(y_rangeNew2);
diff2=length(y_range2)-length(y_rangeNew2)
while (ii~=length(y_range2) && jj<rowy_rangeNew2)
    if (y_rangeNew2(jj)==y_range2(ii))
        jj=jj+1;
        ii=ii+1;
    else
        indexxx=indexxx+1;
        outlier2(indexxx)=ii;
        ii=ii+1;
    end
end
if (indexxx<diff2)
    indexxx=indexxx+1;
    outlier2(indexxx)=length(y_range2)
end
outlierT2 = array2table(outlier2);
%% find the outliers indexs of ymax
jjj=1;
iii=1;
indexxxx=0;
outlier3=[];
rowy_rangeNew3=length(y_rangeNew3);
diff3=length(y_range3)-length(y_rangeNew3);
while (iii~=length(y_range3) && jjj<rowy_rangeNew3)
    if (y_rangeNew3(jjj)==y_range3(iii))
        jjj=jjj+1;
        iii=iii+1;
    else
        indexxxx=indexxxx+1;
        outlier3(indexxxx)=iii;
        iii=iii+1;
    end
end
if (indexxxx<diff3)
    indexxxx=indexxxx+1;
    outlier3(indexxxx)=length(y_range3)
end
outlierT3 = array2table(outlier3);
%% Read data
Y_AfterDelete = table2array(y_rangeT);
Y_Before = y_r(:,6);
Centroid1 = y_r(:,1);
Centroid2 = y_r(:,2);
Major = y_r(:,3);
Orien = y_r(:,5);
%n=length(Y_Before);
r=length(Y_AfterDelete);
w=1;
sumV=0;
%% print all remaining images
i=1;
j=1;
while (i~=r+1&j~=n+1)
    if (Y_AfterDelete(i)==Y_Before(j))
        t = linspace(0,2*pi,50);
        MajorAxisLength(w)=Major(j);
        MinorAxisLength(w)=Y_Before(j);
        a(w) = MajorAxisLength(w)/2;
        b(w) = MinorAxisLength(w)/2;
        Xc(w) = Centroid1(j);
        Yc(w) = Centroid2(j); 
        phi(w) = deg2rad(-Orien(j));
        x = Xc(w) + a(w)*cos(t)*cos(phi(w)) - b(w)*sin(t)*sin(phi(w));
        y = Yc(w) + a(w)*cos(t)*sin(phi(w)) + b(w)*sin(t)*cos(phi(w));
        hold on
        plot(x,y,'k')
        i=i+1; 
        j=j+1;
        w=w+1;
    else
         j=j+1;
    end
end
%% create images for outliers
%find the mean centroid 
centroid1_sum=0;
centroid2_sum=0;
for k=1:n
    centroid1_sum=Centroid1(k)+centroid1_sum;
end
centroid1_mean=centroid1_sum/n;
for k=1:n
    centroid2_sum=Centroid2(k)+centroid2_sum;
end
centroid2_mean=centroid2_sum/n;
%% resize the outliers
lenoutlier=length(outlier);
for k=1:lenoutlier
    if (Y_Before(outlier(k))>30)
        if (outlier(k)==1)
            back=n;
            u=0;
            step=2;
            while (back==outlier(lenoutlier-u))
            back=back-1;
            u=u+1;
            step=step+1;
            end
            front=2;
            uu=2;
            while (front==outlier(uu))
            front=front+1;
            uu=uu+1;
            step=step+1;
            end
            t = linspace(0,2*pi,50);
            distance=abs(Major(front)-Major(back))/step;
            num=outlier(k)-front;
            MajorAxisLength(w)=distance*num+min(Major(front),Major(back));
            a(w) = MajorAxisLength(w)/2;
            MinorAxisLength(w)=Y_Before(front);
            b(w) = MinorAxisLength(w)/2;
            Xc(w) = centroid1_mean;
            Yc(w) = centroid2_mean;
            x = Xc(w) + a(w)*cos(t)*cos(0) - b(w)*sin(t)*sin(0);
            y = Yc(w) + a(w)*cos(t)*sin(0) + b(w)*sin(t)*cos(0);
            hold on
            plot(x,y,'r');
            w=w+1;
            
        elseif(outlier(k)==n)
            back=1;
            u=0;
            step=2;
            while (back==outlier(lenoutlier-u))
            back=back+1;
            u=u+1;
            step=step+1;
            end
            front=n-1;
            uu=2;
            while (front==outlier(uu))
            front=front-1;
            uu=uu+1;
            step=step+1;
            end
            t = linspace(0,2*pi,50);
            distance=abs(Major(front)-Major(back))/step;
            num=outlier(k)-front;
            MajorAxisLength(w)=distance*num+min(Major(front),Major(back));
            a(w) = MajorAxisLength(w)/2;
            MinorAxisLength(w)=Y_Before(front);
            b(w) = MinorAxisLength(w)/2;
            Xc(w) = centroid1_mean;
            Yc(w) = centroid2_mean;
            x = Xc(w) + a(w)*cos(t)*cos(0) - b(w)*sin(t)*sin(0);
            y = Yc(w) + a(w)*cos(t)*sin(0) + b(w)*sin(t)*cos(0);
            hold on
            plot(x,y,'r');
            w=w+1;
        else
            back=outlier(k)+1;
            u=1;
            step=2;
            while ((k+u)<lenoutlier && back==outlier(k+u))
                back=back+1;
                u=u+1;
                step=step+1;
            end
            front=outlier(k)-1;
            uu=1;
            while ( (k-uu>0) && front==outlier(k-uu) )
                front=front-1;
                uu=uu+1;
                step=step+1;
            end
            t = linspace(0,2*pi,40);
            MinorAxisLength(w)=Y_Before(front);
            distance=abs(Major(front)-Major(back))/step;
            num=outlier(k)-front;
            MajorAxisLength(w)=distance*num+min(Major(front),Major(back));
            a(w) = MajorAxisLength(w)/2;
            b(w) = MinorAxisLength(w)/2;
            % Xc = Centroid1(front);
            % Yc = Centroid2(front);
            Xc(w) = centroid1_mean;
            Yc(w) = centroid2_mean;
            x = Xc(w) + a(w)*cos(t)*cos(0) - b(w)*sin(t)*sin(0);
            y = Yc(w) + a(w)*cos(t)*sin(0) + b(w)*sin(t)*cos(0);
            hold on
            plot(x,y,'r');
            w=w+1;
        end
    end
end
%% resize the outliers2
lenoutlier2=length(outlier2);
for k=1:lenoutlier2
    if (Y_Before(outlier2(k))>30)
        if (outlier2(k)==1)
            back=n;
            u=0;
            step=2;
            while (back==outlier2(lenoutlier2-u))
            back=back-1;
            u=u+1;
            step=step+1;
            end
            front=2;
            uu=2;
            while (front==outlier2(uu))
            front=front+1;
            uu=uu+1;
            step=step+1;
            end
            t = linspace(0,2*pi,50);
            distance=abs(Major(front)-Major(back))/step;
            num=outlier2(k)-front;
            MajorAxisLength(w)=distance*num+min(Major(front),Major(back));
            a(w) = MajorAxisLength(w)/2;
            MinorAxisLength(w)=Y_Before(front);
            b(w) = MinorAxisLength(w)/2;
            Xc(w) = centroid1_mean;
            Yc(w) = centroid2_mean;
            x = Xc(w) + a(w)*cos(t)*cos(0) - b(w)*sin(t)*sin(0);
            y = Yc(w) + a(w)*cos(t)*sin(0) + b(w)*sin(t)*cos(0);
            hold on
            plot(x,y,'r');
            w=w+1;
            
        elseif(outlier2(k)==n)
            back=1;
            u=0;
            step=2;
            while (back==outlier2(lenoutlier2-u))
            back=back+1;
            u=u+1;
            step=step+1;
            end
            front=n-1;
            uu=2;
            while (front==outlier2(uu))
            front=front-1;
            uu=uu+1;
            step=step+1;
            end
            t = linspace(0,2*pi,50);
            distance=abs(Major(front)-Major(back))/step;
            num=outlier2(k)-front;
            MajorAxisLength(w)=distance*num+min(Major(front),Major(back));
            a(w) = MajorAxisLength(w)/2;
            MinorAxisLength(w)=Y_Before(front);
            b(w) = MinorAxisLength(w)/2;
            Xc(w) = centroid1_mean;
            Yc(w) = centroid2_mean;
            x = Xc(w) + a(w)*cos(t)*cos(0) - b(w)*sin(t)*sin(0);
            y = Yc(w) + a(w)*cos(t)*sin(0) + b(w)*sin(t)*cos(0);
            hold on
            plot(x,y,'r');
            w=w+1;
        else
            back=outlier2(k)+1;
            u=1;
            step=2;
            while ((k+u)<lenoutlier2 && back==outlier2(k+u))
                back=back+1;
                u=u+1;
                step=step+1;
            end
            front=outlier2(k)-1;
            uu=1;
            while ( (k-uu>0) && front==outlier2(k-uu) )
                front=front-1;
                uu=uu+1;
                step=step+1;
            end
            t = linspace(0,2*pi,40);
            MinorAxisLength(w)=Y_Before(front);
            distance=abs(Major(front)-Major(back))/step;
            num=outlier2(k)-front;
            MajorAxisLength(w)=distance*num+min(Major(front),Major(back));
            a(w) = MajorAxisLength(w)/2;
            b(w) = MinorAxisLength(w)/2;
            % Xc = Centroid1(front);
            % Yc = Centroid2(front);
            Xc(w) = centroid1_mean;
            Yc(w) = centroid2_mean;
            x = Xc(w) + a(w)*cos(t)*cos(0) - b(w)*sin(t)*sin(0);
            y = Yc(w) + a(w)*cos(t)*sin(0) + b(w)*sin(t)*cos(0);
            hold on
            plot(x,y,'r');
            w=w+1;
        end
    end
end
%% resize the outliers3
lenoutlier3=length(outlier3);
for k=1:lenoutlier3
    if (Y_Before(outlier3(k))>30)
        if (outlier3(k)==1)
            back=n;
            u=0;
            step=2;
            while (back==outlier3(lenoutlier3-u))
            back=back-1;
            u=u+1;
            step=step+1;
            end
            front=2;
            uu=2;
            while (front==outlier3(uu))
            front=front+1;
            uu=uu+1;
            step=step+1;
            end
            t = linspace(0,2*pi,50);
            distance=abs(Major(front)-Major(back))/step;
            num=outlier3(k)-front;
            MajorAxisLength(w)=distance*num+min(Major(front),Major(back));
            a(w) = MajorAxisLength(w)/2;
            MinorAxisLength(w)=Y_Before(front);
            b(w) = MinorAxisLength(w)/2;
            Xc(w) = centroid1_mean;
            Yc(w) = centroid2_mean;
            x = Xc(w) + a(w)*cos(t)*cos(0) - b(w)*sin(t)*sin(0);
            y = Yc(w) + a(w)*cos(t)*sin(0) + b(w)*sin(t)*cos(0);
            hold on
            plot(x,y,'r');
            w=w+1;
            
        elseif(outlier3(k)==n)
            back=1;
            u=0;
            step=2;
            while (back==outlier3(lenoutlier3-u))
            back=back+1;
            u=u+1;
            step=step+1;
            end
            front=n-1;
            uu=2;
            while (front==outlier3(uu))
            front=front-1;
            uu=uu+1;
            step=step+1;
            end
            t = linspace(0,2*pi,50);
            distance=abs(Major(front)-Major(back))/step;
            num=outlier3(k)-front;
            MajorAxisLength(w)=distance*num+min(Major(front),Major(back));
            a(w) = MajorAxisLength(w)/2;
            MinorAxisLength(w)=Y_Before(front);
            b(w) = MinorAxisLength(w)/2;
            Xc(w) = centroid1_mean;
            Yc(w) = centroid2_mean;
            x = Xc(w) + a(w)*cos(t)*cos(0) - b(w)*sin(t)*sin(0);
            y = Yc(w) + a(w)*cos(t)*sin(0) + b(w)*sin(t)*cos(0);
            hold on
            plot(x,y,'r');
            w=w+1;
        else
            back=outlier3(k)+1;
            u=1;
            step=2;
            while ((k+u)<lenoutlier3 && back==outlier3(k+u))
                back=back+1;
                u=u+1;
                step=step+1;
            end
            front=outlier3(k)-1;
            uu=1;
            while ( (k-uu>0) && front==outlier3(k-uu) )
                front=front-1;
                uu=uu+1;
                step=step+1;
            end
            t = linspace(0,2*pi,40);
            MinorAxisLength(w)=Y_Before(front);
            distance=abs(Major(front)-Major(back))/step;
            num=outlier3(k)-front;
            MajorAxisLength(w)=distance*num+min(Major(front),Major(back));
            a(w) = MajorAxisLength(w)/2;
            b(w) = MinorAxisLength(w)/2;
            % Xc = Centroid1(front);
            % Yc = Centroid2(front);
            Xc(w) = centroid1_mean;
            Yc(w) = centroid2_mean;
            x = Xc(w) + a(w)*cos(t)*cos(0) - b(w)*sin(t)*sin(0);
            y = Yc(w) + a(w)*cos(t)*sin(0) + b(w)*sin(t)*cos(0);
            hold on
            plot(x,y,'r');
            w=w+1;
        end
    end
end
%% Calculate Error
for q=1:(w-1)
    v(q)=4/3*pi * a(q) * b(q) * b(q);
end

vol=reshape(v,w-1,1);
vv=1;
for a_=40:70
    x=3.7/a_;
    c=1;
    SumV=0;
    while (c~=w)
        volume(c,vv)=vol(c) * x^3;
        SumV=SumV+volume(c,vv);
        c=c+1;
    end
    ave(index,vv)=SumV/(w-1);
    error(index,vv)=abs( ave(index,vv) - 90)/90 * 100;
    vv=vv+1;
end 
  