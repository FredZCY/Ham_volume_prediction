%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This script is to perform the defect detection on ham images
%
%  Reference:
%  xxxxx
%
%  Matlab version was written and was tested on Matlab 2020a
%  If you have any problem, please feel free to contact Sze Teng Liong (stliong@fcu.edu.tw)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;close all;clear
%%
index=1
str=dir(['C:\Users\lenovo\Desktop\HAM\Hum5\test_v_maskoutput',num2str(index,'%02d'),'b\*.jpg'])
[n,~]=size(str)
for i=1:n
    m1=imread(['C:\Users\lenovo\Desktop\HAM\Hum5\test_v_maskoutput',num2str(index,'%02d'),'b\p',num2str(i),'.jpg']); 
    BW1 = im2bw(m1);
    se = strel('disk',100);
    BW2 = imclose(BW1,se);
    s{i} = regionprops('table',BW2,{...
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
minorr=minor_new(:,4);
%% delte extremely small ones
numRow=length(minorr);

m=max(minorr);
numSmall=numRow-n;
minorr_New=minorr;

%% Find out the IQR
Q1=quantile(minorr_New,0.25,1);
Q3=quantile(minorr_New,0.75,1);
IQR = iqr(minorr_New);
high=Q3+1.5*IQR;
low=Q3-1.5*IQR;
bound_up=Q3+3*IQR;
bound_low=Q1-3*IQR;
%% Delete those exceed 1.5IQR     
outlierType1=minorr_New>high|minorr_New<low;
outlierType2=minorr_New>bound_up|minorr_New<bound_low;
minorr_New(outlierType1)=[];
while (any(outlierType2==1))
    Q1=quantile(minorr_New,0.25,1);
    Q3=quantile(minorr_New,0.75,1);
    IQR = iqr(minorr_New);
    top=Q3+1.5*IQR;
    low=Q3-1.5*IQR;
    bound_up=Q3+3*IQR;
    bound_low=Q1-3*IQR;
    outlierType1=minorr_New>top|minorr_New<low;
    outlierType2=minorr_New>bound_up|minorr_New<bound_low;
    minorr_New(outlierType1)=[];
end
%% save MinorAxisLength after delete all outliers
minor_T= array2table(minorr_New);

%% find the outliers indexs
j=1;
i=1;
indexx=0;
row_minor_New=length(minorr_New);
diff=length(minorr)-length(minorr_New)
while (i~=length(minorr) || j<row_minor_New)
    if (minorr_New(j)==minorr(i))
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
    outlier(indexx)=length(minorr)
end
outlierT = array2table(outlier);

%% Read data
Minor_AfterDelete = table2array(minor_T);
Minor_Before = minor_new(:,4);
Centroid1 = minor_new(:,1);
Centroid2 = minor_new(:,2);
Major = minor_new(:,3);
Orien = minor_new(:,5);
n=length(Minor_Before);
r=length(Minor_AfterDelete);
w=1;
sumV=0;
%% print all remaining images
i=1;
j=1;
while (i~=r+1&j~=n+1)
    if (Minor_AfterDelete(i)==Minor_Before(j))
        t = linspace(0,2*pi,50);
        MajorAxisLength(w)=Major(j);
        MinorAxisLength(w)=Minor_Before(j);
        a(w) = MajorAxisLength(w)/2;
        b(w) = MinorAxisLength(w)/2;
        Xc(w) = Centroid1(j);
        Yc(w) = Centroid2(j); 
        phi(w) = deg2rad(-Orien(j));
        x = Xc(w) + a(w)*cos(t)*cos(phi(w)) - b(w)*sin(t)*sin(phi(w));
        y = Yc(w) + a(w)*cos(t)*sin(phi(w)) + b(w)*sin(t)*cos(phi(w));
        hold on
        plot(x,y,'k');
        i=i+1; 
        j=j+1;
        w=w+1;
    else
         j=j+1;
    end
end
%% create images for outliers
% outlier = array2table(outlierr);
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
    if (Minor_Before(outlier(k))>30)
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
            MinorAxisLength(w)=Minor_Before(front);
            a(w) = MajorAxisLength(w)/2;
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
            MinorAxisLength(w)=Minor_Before(front);
            a(w) = MajorAxisLength(w)/2;
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
            MinorAxisLength(w)=Minor_Before(front);
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
