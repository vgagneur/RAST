function [coefmat,cuth,cutl]=RAST_picture_segmentation(I,ii,no,path,pixsize,numlines,degrot)

%Original version by Vincent Gagneur.
%This function segments the picture in small squares of pixsize*pixsize px dimension.
%it then performs an edge detection in each segment.
%it then finds n=numlines lines in each segment using a hough transform
%algorithm.
%The line detected (slopes and coordinates) are stored in "coefmat" a
%matrix for (hopefully) fast manipulation.

[h,l]=size(I);
seg=0;
cuth=0;


%finds how many segments can fit in the picture, could be made faster but it's already quite fast anyway
for i=1:1:h
    seg=seg+pixsize;
    if seg<h
       cuth=cuth+1  ;
        %pixsizexpixsize pixel for segmented areas
    else
        break
    end
end

seg=0;
cutl=0;

for i=1:1:l
    seg=seg+pixsize;
    if seg<l
       cutl=cutl+1  ;
        %pixsizexpixsize pixel for segmented areas
    else
        break
    end
end

%[h2,l2]=size(I); %just for a check
%adding a black dot to the corners to avoid them. (prob obsolete)
I(:,1)=0;
I(1,:)=0;
I(:,l)=0;
I(h,:)=0;

while 1

a=1;
b=1;

c=pixsize; %idivide(int16(h),int16(cut));
d=pixsize; %idivide(int16(l),int16(cut));

e=pixsize; %idivide(int16(h),int16(cut));
f=pixsize; %idivide(int16(l),int16(cut));

angleList = -90:0.1:89; %number of possible slopes
%angleList(angleList==45 | angleList==-45 | angleList==-90 | angleList==90) = [];% trying to get rid of some mathematical errors by removing 45 and 90° lines (no longer effective when the angle step is very small)

%cut loop:

L=I;

imname=sprintf('%.0f_picture.png', ii);
locname=fullfile(path, imname);
imwrite(I, locname);

Li=L;
J=L(a:c,b:d);

cutcount=0;
cutmax=cuth*cutl;

%g = waitbar(0,sprintf('Segmenting picture and searching for lines (picture %.0f / %.0f)',[ii no])); %note that waitbars slow down the code quite a lot, so try and remove them if you care about speed
fprintf('Segmenting picture and searching for lines (picture %.0f / %.0f)',[ii no])
fprintf('\n');

for i = 1:1:(cuth)
    for j = 1:1:(cutl)

        cutcount=cutcount+1;
        %waitbar(cutcount/cutmax)

        J=L(a:c,b:d);
        [xj,yj]=size(J);
        Ji=J(2:xj-1,2:yj-1); %getting rid of the edges of the segment
% 
        Ji=edge(Ji,'Canny'); % Edge detection step for figure making % note: obsolete??

        if ismember(0,J)==0 %will ignore area if part of masked areas, which means 0 is present
            
            J=imadjust(J); %auto contrast using imadjust, honestly not quite sure if it matters

            [xj,yj]=size(J);
            Ji=J(2:xj-1,2:yj-1); %getting rid of the edges
            
            Ji=edge(Ji,'Canny'); %edge locally for line detection

            %detecting lines in each segments using hough transform
            [H,theta,rho] = hough(Ji, 'Theta', angleList);
            P = houghpeaks(H,numlines,'threshold',0.3*max(H(:)));
            lines = houghlines(Ji,theta,rho,P,'FillGap',100,'MinLength',4); 
            %houghlines returns a structure, the next bit expands it to
            %store all lines detected in each segments.
            %initialisiation of line storage:
            xycell={};
            wzcell={};
            if exist('alllines','var')==0 %initialize line storage
                for ich=1:1:length(lines)
                    xycell{ich}=double([b,a]);
                    wzcell{ich}=double([d,c]);
                    icell{ich}=i;
                    jcell{ich}=j;
                end

                alllines=lines;
                [alllines.xyco]=xycell{:};
                [alllines.wzco]=wzcell{:};
                [alllines.celx]=icell{:};
                [alllines.cely]=jcell{:};

            else
            %else stores directly
                for ll = 1:1:length(lines) %if chooses to find multiple line per area
                    for ich=1:1:length(lines)
                        xycell{ich}=double([b,a]);
                        wzcell{ich}=double([d,c]);
                        icell{ich}=i;
                        jcell{ich}=j;
%                         neibcell{ich}=0;
                    end
                    [lines.xyco]=xycell{:};
                    [lines.wzco]=wzcell{:};
                    [lines.celx]=icell{:};
                    [lines.cely]=jcell{:};
%                     [lines.neib]=neibcell{:};
%                     [lines.best]=neibcell{:};
%                     [lines.second]=neibcell{:};
%                     [lines.coef]=neibcell{:};
                    alllines(end+1)=lines(ll);
                    
                end
            end
        end

        L(a+1:c-1,b+1:d-1)=Ji;
        Li(a:c,b:d)=J;
        %imshow(Li)
        b=b+f;
        d=d+f;
    end
    
    b=1;
    %b=round(h/cut)-1;
    df=d-f;
    d=pixsize; %idivide(int16(l),int16(cut));
    cf=c;
    a=a+e;
    c=c+e;
end

%close(g)
pause(0.1)

a=1;
b=1;
L=L(a:cf,b:df); %taking the "final coordinates" of the picture
L=(L>0);
BW=L;
    for k = 1:length(alllines)
       %replacing points coordinates by new coordinates

       alllines(k).point1=alllines(k).point1+alllines(k).xyco;
       alllines(k).point2=alllines(k).point2+alllines(k).xyco;

    end
    
    break 
end
clf

imname=sprintf('%.0f_contrast.png', ii);
locname=fullfile(path, imname);
imwrite(Li, locname); %saves contrast adjusted picture

imname=sprintf('%.0f_bw.png', ii);
locname=fullfile(path, imname);
imwrite(BW, locname); % saves edge detection figure

hold off %useless?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates a matrix of coefficients, coordinates for later use. Matrix is
% usefull, it allows quicker computing (maybe)

coefmat=zeros(cuth,cutl,7*numlines);
k=0;

%for some reason the first segment has 5 lines attributed to it, so here it removes 2
%of them (this shows the risk of this storing technique, using cells may
%be more reliable. update pre-publication: never had an issue in 3 years of
%use.

alllines=alllines(3:end);

while k <=(length(alllines)-numlines)
    for l=1:1:numlines
       k=k+1;
       %xy = [alllines(k).point1; alllines(k).point2];
       A=alllines(k).point1;
       B=alllines(k).point2;
       cur_deg=atand((A(2)-B(2))/(B(1)-A(1))); %slope in degrees
       X=alllines(k).celx;
       Y=alllines(k).cely;
       coefmat(X,Y,l)=cur_deg+degrot; % Stores line orientation. If degrot =/= 0, rotates the lines if EBSD map at different rotation than SEM picture, theoretical slip plots are also rotated by degrot on the figures.
       coefmat(X,Y,numlines+1+6*(l-1))=A(1); %Ax
       coefmat(X,Y,numlines+2+6*(l-1))=A(2); %Ay
       coefmat(X,Y,numlines+3+6*(l-1))=B(1); %Bx
       coefmat(X,Y,numlines+4+6*(l-1))=B(2); %By
       coefmat(X,Y,numlines+5+6*(l-1))=0; %will store line score 
       coefmat(X,Y,numlines+6+6*(l-1))=0; %will store best match for figure making

    end
end

%coefmat;

close all

end