%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
%
% RAST: Rapid Analysis of Slip Traces, original version written by Vincent Gagneur (2022-2025)
% The release of the code accompanies a paper which describes its main
% functionalities, and how to interpet the data output in more details, as applied to indentation generated slip traces for bcc/B2 slip activities,  see V. Gagneur et. al. "automated indentation-based slip trace
% analysis for bcc and B2 plasticity" (expected publication in 2025)
%
% This code requires the "image processing toolbox" from mathworks (install
% via "apps", "get more apps", and find it there). It is possible that it
% depends on other functions from other apps (which I forgot about), just read the error messages when running the
% code and it should be indicated which matlab toolbox you need to install.
%
% This matlab file is used to load a new dataset, prompting the user to
% define each picture, crop/mask unwanted features, and attribute Euler
% angles to each picture.
% 
% In this part of the code you have:
%
% (1) picture loading by user,
%
% (2)/(3) picture cropping / masking by user, to hide indent, other grain
% areas, and large highly contrasted linear feature that risk getting detected
%
% (4) User manually associating the EBSD determined Euler angles.
% 
% Important!!!, currently pictures from the same grain orientation should
% be entered one after the others, or the grain ID numeration in "(8)"
% output by the code will mistakenly give different grain numbers to pictures from the same
% grain. (note that the actual Euler angle values are output in the excel
% file, so user can double check here whether the grain IDs are correctly attributed and manually fix it if needed).
%
% 
%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

prompt = {'Chose a name for your file. It will create a folder where pictures and results are stored'}; %just giving a name to your folder where everything (figures and excel spredsheet) will be stored) It also give the name to the .mat file created for re-use of pictures and Eulers inputted. Matlab usually doesnt like spaces between characters
dlgtitle = 'Saving the table';
dims = [1 100];
filename = inputdlg(prompt,dlgtitle,dims);
%mkdir(filename{1})
path=strcat(pwd,'\',filename{1});

no=0;
Eulers={};
Pictures={};
BWss={};
xi2={};
yi2={};

%%%%%%%%%%%%%%%%%%%%
% (1) Picture loading and indent masking loop

while 1
    
    no=no+1;
    
    %Picture choice + crop + mask loop:

    while 1 

        file = uigetfile('*.*','Select your picture');
        I = imread(file); % <---- your picture must be in the same folder as the code
        I = im2gray(I);

        f=figure; 
        hold on
        imshow(I)

        %%%%%%%%
        % (2) cropping: (optional, used mainly to crop out SEM image legend (bottom e.g. magnification, voltage etc.))

        imshow(I)

        dim = [0.2 0.6 .3 .3];
        str = 'Crop out unwanted parts of the picture, drag area to keep, right click crop when happy';
        a=annotation('textbox',dim,'String',str,'FitBoxToText','on');
        a.FontSize = 12;
        a.Color = 'r';
        a.FontWeight = 'bold';
        a.LineStyle='none';

        I=imcrop; 

        %%%%%%%%%
        % (3) Manual Masking:

        imshow(I)

        dim = [0.2 0.6 .3 .3];
        str = 'Draw polygon to mask indent: Left click to draw segments, Right click to complete, then right click and select create mask at the end'; %note to fix, text doesn't show up i believe
        aa=annotation('textbox',dim,'String',str,'FitBoxToText','on');
        aa.FontSize = 12;
        aa.Color = 'r';
        aa.FontWeight = 'bold';
        aa.LineStyle='none';

        [xrefout,yrefout,BWs,xi,yi]= roipoly(I); %masking indent

        I(BWs)=0.5; %put at 0.5 for the smoothing

        imshow(I)
        
        continu=questdlg('Do you need to re-crop, re-mask or use a different picture? (this picture will be forgoten!)','continue prompt') ;
        if isequal(continu,'No')
           break;
        end  

    end

    %%%%%%%%%%%%%
    % (4) Euler prompt:
    prompt = {'Euler format [x,y,z]:'}; % e.g. type in the prompt: [90,50,131]
    dlgtitle = 'Enter the Euler for this picture';
    dims = [1 100];
    Euler = inputdlg(prompt,dlgtitle,dims);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Storing both the euler and the picture, and the masks:
    Eulers{no}=Euler;
    Pictures{no}=I;
    BWss{no}=BWs;
    xi2{no}=xi;
    yi2{no}=yi;
    
    continu=questdlg('Do you have other orientations/pictures to enter for this sample?','continue prompt') ;
    if isequal(continu,'No')
        break;
    end  

end

close all

save(filename{1},"no","Pictures","Eulers","BWss","xi2","yi2") %saving all above as a .mat file (same folder as code) for quick re-use. To perform the analysis, use the final name as input in RAST_Running.m 
