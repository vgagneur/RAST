function  []=RAST_Figuremaking(I,ii,nei,pixsize,error,cuth,cutl,numlines,coefmat,Euler,degrot,path,crystal)
%%%%%%%%%%%%%%
% The use of this code requires another code "Euler_to_Slip.m" which
% calculates the theoretical slip traces from Eulers angles.
%
% This code simply plots the theoretical orientations for each plane, for
% figure making, degrot accounts for user inputed rotation between SEm and EBSD imaging if known.
%
%%%%%%%%%%%%%%

figfin=figure('units','normalized','outerposition',[0 0 0.8 1]);
pause(0.1)
tiledlayout(6,9)
nexttile(1,[6 6])
imshow(I);
axis on

[t,s]=title(sprintf('Picture %.0f: lines with a score >= %.0f',[ii nei]),...
    ['euler = ',Euler{1},' ','/',' N = ',num2str(pixsize),'*',num2str(pixsize),' px ','/',' IR = ', num2str(error),'°', ' rotation correction = ', num2str(degrot),'°']);
s.FontAngle = 'italic';
axis off
hold on

for i = 1:1:(cuth)
    for j = 1:1:(cutl)
        for k=1:1:numlines
            if coefmat(i,j,numlines+5+6*(k-1))>=nei 

                A=[coefmat(i,j,numlines+1+6*(k-1)),coefmat(i,j,numlines+2+6*(k-1))]; 
                B=[coefmat(i,j,numlines+3+6*(k-1)),coefmat(i,j,numlines+4+6*(k-1))];

                xy = [A; B];

                %colors the line depending on its best matching
                %theoretical plane, or magenta if unambiguously
                %<111> (did not match 100 nor 110)

                coefmat(i,j,numlines+6+6*(k-1));

                if  coefmat(i,j,numlines+6+6*(k-1))==100
                    plot(xy(:,1),xy(:,2),'-','LineWidth',1,'Color',[0,0,255/255]); % (100) in blue

                elseif coefmat(i,j,numlines+6+6*(k-1))==110
                    plot(xy(:,1),xy(:,2),'-','LineWidth',1,'Color',[255/255,128/255,0]); % (110) in orange

                elseif coefmat(i,j,numlines+6+6*(k-1))==112
                    plot(xy(:,1),xy(:,2),'-','LineWidth',1,'Color',[0,255/255,0]); % (112) in green

                elseif coefmat(i,j,numlines+6+6*(k-1))==123
                    plot(xy(:,1),xy(:,2),'-','LineWidth',1,'Color',[255/255,255/255,0]); % (123) in yellow

                elseif coefmat(i,j,numlines+6+6*(k-1))==111 % special condition in magenta, defined in line_processing.m, for bcc unambiguously <111>,  for fcc ambiguous match
                    plot(xy(:,1),xy(:,2),'-','LineWidth',1,'Color','m');

                elseif coefmat(i,j,numlines+6+6*(k-1))==0 % if matches nothing, show in red
                    plot(xy(:,1),xy(:,2),'-','LineWidth',1,'Color','red'); 

                end

            end

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uses the Euler value and calculates intersection using Euler_to_Slip.m

[slip100,slip110,slip112,slip123]=RAST_euler_to_slip(Euler,crystal); % <-- Euler_to_Slip.m must be in the same folder.
%normalisation of vectors
slips = {slip100, slip110, slip112, slip123};

for j=1:length(slips)
    slipette=slips{j}; 
    for i=1:length(slipette)
        slipette(i,:)=slipette(i,:)/norm(slipette(i,:));
    end
    slips{j}=slipette;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The next loop converts the vectors (slip directions) to corresponding
%slopes

CO={}; %matrix where corresponding slopes will be stored
for i = 1:length(slips)
    slip=transpose(slips{i});
    coefstore=zeros(1,length(slip));
    for j = 1:length(slip)
        slope=(slip(2,j))/(-slip(1,j)); %determines the slope of each slip, slope = y / x, correction x --> -x 
        coefstore(j)=slope;
    end
    CO{i}=coefstore;
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plots pole figures with expected slip trace orientations

for i = 1:length(slips)
    slip=transpose(slips{i});
    hold on
    nexttile(9*(i+1)-2)
    %subplot(2,2,i);

    if crystal==1
        names={'{100}','{110}','{112}','{123}'};
    elseif crystal==2
        names={'{111}','{-111}','{1-11}','{11-1}'};
    end
    
    
    for j = 1:length(slip) %( -slip x)

        a=(slip(1,j));
        b=(slip(2,j));

        if degrot==0
            x=a;
            y=b;
            %nothing
        else %rotate vectors
            x=a*cosd(degrot)-b*sind(degrot); % NB +degrot works (should have been -degrot), probably because slip coordinates are inverted compared to the picture axis so gets inverted when plotted on graph / figures.
            y=a*sind(degrot)+b*cosd(degrot);
        end

        if i==1
            quiver(0,0,-x,y,'Color',[0,0,255/255]);
        elseif i==2
            quiver(0,0,-x,y,'Color',[255/255,128/255,0]);
        elseif i==3
            quiver(0,0,-x,y,'Color',[0/255,255/255,0/255]);
        elseif i==4
            quiver(0,0,-x,y,'Color',[230/255,230/255,0]);
        end

        hold on

    end

    %plots circles and titles
    title(names{i});

    if i==1
        viscircles([0,0],1,'Color',[0,0,255/255]);
    elseif i==2
        viscircles([0,0],1,'Color',[255/255,128/255,0]);
    elseif i==3
        viscircles([0,0],1,'Color',[0/255,255/255,0/255]);
    elseif i==4
        viscircles([0,0],1,'Color',[250/255,250/255,0]);
    end
    
    axis equal
    xlim([-1 1])
    ylim([-1 1])
    axis off
    a=1;

end

hold off

%%%%%%%%%%%%
%Saving figure

imname=sprintf('picture_%.0f_ST_%.0f_linesselected.png', [ii nei]);
locname=fullfile(path, imname);
exportgraphics(figfin, locname,'Resolution',200)
%saveas(figfin, locname);

hold off

close all
pause(0.1)



    
end
    


            
            
            
            
            
            


            
    