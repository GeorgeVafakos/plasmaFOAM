function [ ] = cylindricalMeshGenerator( )
% Creates structured hex butterfly mesh for OpenFOAM blockMeshDict

% manual uses the command window to enter input data, auto uses the input data in the .m file
    rData= input('0 - manual input, 1 - auto input: ');

    if rData
        % Give array of coordinates along the symmetry axis (in [ ... ] format)
        z = [-0.1025,-0.0025,0.0025,0.1025];
        % Give array of coordinates along the radial axis, excluding r=0 (in [ ... ] format)
        r = [0.0008,0.0016,0.0032,0.0034,0.010];
        % Give array of BLOCKS to be excluded (in [ ... ] format)
        B = [10];
        % Give axial mesh distribution
        MSHz = [80,80,80];
        % Give axial mesh grading
        GRDNGz = [0.01,1,100];
        % Give radial mesh distribution
        MSHr = [10,10,10,5,20];
        % Give radial mesh grading
        GRDNGr = [1,1,1,1,20];
        % Give convertToMeters
        convertToMeters = 1;
        % Give arc rate of 1st block, default is 0.5 (1 for circle, 0 for square)
        rate1st = 0.5;
        % Complete the 'Boundary' section in the blockMeshDict?
        iRunBs = 1;
        % Define up to   9   Boundaries
        % Give names of Boundaries
        name1 = 'leftPlateHe';
        name2 = 'leftPlateDiel';
        name3 = 'leftPlateAir';
        name4 = 'airOuter';
        name5 = 'electrode';
        name6 = 'rightPlateHe';
        name7 = 'rightPlateDiel';
        name8 = 'rightPlateAir';
        % Give types of Boundaries, 1 - wall, 2 - inlet, 3 - outlet, 4 - patch
        types = [4,4,4,4,4,4,4,4];
        % Give faces of Boundaries
        bounds1 = [1,2];
        bounds2 = [3];
        bounds3 = [4,6];
        bounds4 = [13,14,16];
        bounds5 = [5,10,15,17];
        bounds6 = [7,8];
        bounds7 = [9];
        bounds8 = [11,12];


        iz=length(z);
        ir=length(r);
    else
        rData=0;
        z=input('Give array of coordinates along the symmetry axis (in [ ... ] format): ');
        r=input('Give array of coordinates along the radial axis, excluding r=0 (in [ ... ] format): ');

        iz=length(z);
        ir=length(r);

    %Geometry
        plotBlocks(z,r);
        B=input('Give array of BLOCKS to be excluded (in [ ... ] format): ');



    %Meshing

        plotMesh(z,r);
        [ MSHz , GRDNGz , MSHr , GRDNGr ] = createMesh(z,r);


        convertToMeters = input('Give convertToMeters, default is 1: ');
        if isempty( convertToMeters )
            convertToMeters = 1;
        end
    end

    fID=fopen('blockMeshDict','w');
    fprintf(fID,'/*--------------------------------*- C++ -*----------------------------------*\\ \n');
    fprintf(fID,'  =========                 |\n');
    fprintf(fID,'  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n');
    fprintf(fID,'   \\\\    /   O peration     | Website:  https://openfoam.org\n');
    fprintf(fID,'    \\\\  /    A nd           | Version:  7\n');
    fprintf(fID,'     \\\\/     M anipulation  |\n');
    fprintf(fID,'\\*---------------------------------------------------------------------------*/\n');
    fprintf(fID,'FoamFile\n{\nversion     2.0;\nformat      ascii;\nclass       dictionary;\nobject      blockMeshDict;\n}\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\nconvertToMeters %8.4e;\n\n',convertToMeters);

    % Points

    fprintf(fID,'vertices\n(\n');
    count=0;
    s2=sqrt(2);
    for i=1:ir
        for j=1:iz
            fprintf(fID,'(%12.8e  %12.8e %12.8e)  // p%3i \n',-r(i)/s2,-r(i)/s2,z(j),count);
            count=count+1;
            fprintf(fID,'(%12.8e  %12.8e %12.8e)  // p%3i \n', r(i)/s2,-r(i)/s2,z(j),count);
            count=count+1;
            fprintf(fID,'(%12.8e  %12.8e %12.8e)  // p%3i \n', r(i)/s2, r(i)/s2,z(j),count);
            count=count+1;
            fprintf(fID,'(%12.8e  %12.8e %12.8e)  // p%3i \n',-r(i)/s2, r(i)/s2,z(j),count);
            count=count+1;
        end
    end

    %Blocks
    fprintf(fID,');\n\nblocks\n(\n');
    count=0;
    for j=1:iz-1
        if isempty(find(B==count))
            base = count;
            fprintf(fID,'    hex (%3i  %3i  %3i  %3i  %3i  %3i  %3i  %3i)( %4i %4i %4i ) simpleGrading ( %2i %2i %2i ) // block %2i \n',4*base,4*base+1,4*base+2,4*base+3,4*base+4,4*base+5,4*base+6,4*base+7,2*MSHr(1),2*MSHr(1),MSHz(j),1,1,GRDNGz(j),count);
        end
        count=count+1;
    end
    count1_4=count;
    for i=1:ir-1
        for j=1:iz-1
            if isempty(find(B==count1_4))
                base1 = iz*(i-1)*4+(j-1)*4;
                base2 = iz*i*4+(j-1)*4;
                fprintf(fID,'    hex (%3i  %3i  %3i  %3i  %3i  %3i  %3i  %3i)( %4i %4i %4i ) simpleGrading ( %2i %2i %2i ) // block %2i \n',base1+1,base1,  base2,  base2+1,base1+5,base1+4,base2+4,base2+5,2*MSHr(1),MSHr(i+1),MSHz(j),1,GRDNGr(i+1),GRDNGz(j),count);
                count=count+1;
                fprintf(fID,'    hex (%3i  %3i  %3i  %3i  %3i  %3i  %3i  %3i)( %4i %4i %4i ) simpleGrading ( %2i %2i %2i ) // block %2i \n',base1+2,base1+1,base2+1,base2+2,base1+6,base1+5,base2+5,base2+6,2*MSHr(1),MSHr(i+1),MSHz(j),1,GRDNGr(i+1),GRDNGz(j),count);
                count=count+1;
                fprintf(fID,'    hex (%3i  %3i  %3i  %3i  %3i  %3i  %3i  %3i)( %4i %4i %4i ) simpleGrading ( %2i %2i %2i ) // block %2i \n',base1+3,base1+2,base2+2,base2+3,base1+7,base1+6,base2+6,base2+7,2*MSHr(1),MSHr(i+1),MSHz(j),1,GRDNGr(i+1),GRDNGz(j),count);
                count=count+1;
                fprintf(fID,'    hex (%3i  %3i  %3i  %3i  %3i  %3i  %3i  %3i)( %4i %4i %4i ) simpleGrading ( %2i %2i %2i ) // block %2i \n',base1  ,base1+3,base2+3,base2,  base1+4,base1+7,base2+7,base2+4,2*MSHr(1),MSHr(i+1),MSHz(j),1,GRDNGr(i+1),GRDNGz(j),count);
                count=count+1;
            end
            count1_4=count1_4+1;
        end
    end

    %Edges
    fprintf(fID,');\n\nedges\n(\n');
    if ~rData
        rate1st = input('Give arc rate of 1st block, default is 0.5 (1 for circle, 0 for square): ');
    end
    if isempty(rate1st);
        rate1st=0.5;
    end
    for j=1:iz
        base = (j-1)*4;
        arc=r(1)*(rate1st+(1-rate1st)/s2);
        fprintf(fID,'    arc %3i  %3i  (%12.8e  %12.8e %12.8e)\n', base  , base+1,0,-arc,z(j));
        fprintf(fID,'    arc %3i  %3i  (%12.8e  %12.8e %12.8e)\n', base+1, base+2, arc,0,z(j));
        fprintf(fID,'    arc %3i  %3i  (%12.8e  %12.8e %12.8e)\n', base+2, base+3,0, arc,z(j));
        fprintf(fID,'    arc %3i  %3i  (%12.8e  %12.8e %12.8e)\n', base+3, base  ,-arc,0,z(j));
    end

    for i=2:ir
        for j=1:iz
            base = iz*(i-1)*4+(j-1)*4;
            arc=r(i);%r(1)*(rate1st+(1-rate1st)/s2);
            fprintf(fID,'    arc %3i  %3i  (%12.8e  %12.8e %12.8e)\n', base  , base+1,0,-arc,z(j));
            fprintf(fID,'    arc %3i  %3i  (%12.8e  %12.8e %12.8e)\n', base+1, base+2, arc,0,z(j));
            fprintf(fID,'    arc %3i  %3i  (%12.8e  %12.8e %12.8e)\n', base+2, base+3,0, arc,z(j));
            fprintf(fID,'    arc %3i  %3i  (%12.8e  %12.8e %12.8e)\n', base+3, base  ,-arc,0,z(j));
        end
    end

    %Boundary

    if ~rData
        iRunBs=input('0 - stop, 1 - give boundaries: ');
        if isempty(iRunBs)
            iRunBs=0;
        end
    end

    fprintf(fID,');\n\nboundary\n(\n');

    [Bs,nBs]=Boundaries(z,r,B);

    if iRunBs

        checkBs=1:nBs;
        type=['wall  ';'inlet ';'outlet';'patch '];
        count=0;

        while find(checkBs)
            count=count+1;
            if ~rData
                name = input('Give boundary name (in ''...''): ');
                itype = input('Choose BC type 1 - wall, 2 - inlet, 3 - outlet, 4 - patch : ');
                while itype<1 || itype>4
                    itype = input('Incorrect Type!! Choose BC type 1 - wall, 2 - input, 3 - output, 4 - patch: ');
                end
                bounds = input('Give array of boundary faces (in [ ... ] format): ');
            else
                itype=types(count);
                switch count
                    case 1
                        name = name1; bounds = bounds1;
                    case 2
                        name = name2; bounds = bounds2;
                    case 3
                        name = name3; bounds = bounds3;
                    case 4
                        name = name4; bounds = bounds4;
                    case 5
                        name = name5; bounds = bounds5;
                    case 6
                        name = name6; bounds = bounds6;
                    case 7
                        name = name7; bounds = bounds7;
                    case 8
                        name = name8; bounds = bounds8;
                    case 9
                        name = name9; bounds = bounds9;
                end
            end

            for i=1:length(bounds)
                checkBs(find(checkBs == bounds(i)))=0;
            end

            fprintf(fID,'   %s\n   {\n      type %s;\n      faces (\n',name,type(itype,:));
            for i=1:length(bounds)
                if isempty(find(Bs(bounds(i),5:12)))
                    fprintf(fID,'         (%3i %3i %3i %3i)\n',Bs(bounds(i),1:4));
                else
                    fprintf(fID,'         (%3i %3i %3i %3i)\n         (%3i %3i %3i %3i)\n         (%3i %3i %3i %3i)\n         (%3i %3i %3i %3i)\n',Bs(bounds(i),:));
                end
            end
            
            if ~rData
                disp('Remaining boundary faces:')            
                disp(checkBs(find(checkBs)))
            end

            fprintf(fID,'      );\n   }\n');
        end
    end
    fprintf(fID,');\n\nmergePatchPairs();');
    fclose(fID);

end
function [  ] = plotBlocks(z,r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    iz=length(z);
    ir=length(r);
    close all;
    figure;
    axis([min(z)-0.1*abs(max(z)) , max(z)+0.1*abs(max(z)), 0 , max(r)+0.1*abs(max(r))]);
    hold;
    
    x=[min(z),max(z)];
    for i=1:ir
        y=[r(i),r(i)];
        plot(x,y);
    end
    y=[0,max(r)];
    for i=1:iz
        x=[z(i),z(i)];
        plot(x,y);
    end

    count=0;
    for j=1:iz-1
        text((z(j)+z(j+1))/2,(0+r(1))/2,num2str(count),'Color','red');
        count=count+1;
    end
    for i=1:ir-1
        for j=1:iz-1
            text((z(j)+z(j+1))/2,(r(i)+r(i+1))/2,num2str(count),'Color','red');
            count=count+1;
        end
    end
    hold;

end

function [  ] = plotMesh( z,r )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    iz=length(z);
    ir=length(r);
    close all;
    figure;
    axis([min(z)-0.1*abs(max(z)) , max(z)+0.1*abs(max(z)), 0 , max(r)+0.1*abs(max(r))]);
    hold;
    
    x=[min(z),max(z)];
    for i=1:ir
        y=[r(i),r(i)];
        plot(x,y);
    end
    y=[0,max(r)];
    for i=1:iz
        x=[z(i),z(i)];
        plot(x,y);
    end
    count=0;
    for j=1:iz-1
        text((z(j)+z(j+1))/2,0,num2str(count),'Color','red');
        count=count+1;
    end
    text(z(1),r(1)/2,num2str(count),'Color','red');
    count=count+1;
    for i=1:ir-1
        text(z(1),(r(i)+r(i+1))/2,num2str(count),'Color','red');
        count=count+1;
    end
    hold;

end


function [ MSHz , GRDNGz , MSHr , GRDNGr ] = createMesh( z,r )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    iz=length(z);
    ir=length(r);
    
    expon = 0.75;
    iniEst=50;
    estimate = round((z(2)-z(1))/max(z)*iniEst)+1;
    MSHz = input(['Give mesh size for face 0, default is ',num2str(estimate),': ']);
    if isempty(MSHz)
        MSHz=estimate;
    end
    GRDNGz = input(['Give grading for faces 0 - ',num2str(iz-2),',default is [ones] (in [ ... ] format): ']);
    if isempty(GRDNGz)
        GRDNGz = ones(1,iz-1);
    end
    while length(GRDNGz)~=iz-1
        GRDNGz = input(['Icorrect array length!! Give grading for faces 0 - ',num2str(iz-2),',default is [ones] (in [ ... ] format): ']);
    end
    for i=2:iz-1
        MSHz(i)= round(MSHz(i-1)*(z(i+1)-z(i))/(z(i)-z(i-1))/GRDNGz(i)^expon)+1;
    end
    disp('Estimated axial mesh distribution:')
    disp(MSHz)
    temp = input('Update axial mesh distribution, or press ENTER to continue: ');
    if ~isempty(temp)
        MSHz=temp;
    end

    estimate = round(r(1)/max(r)*(max(r)/max(z)*iniEst))+1;
    MSHr = input(['Give mesh size for face 0, default is ',num2str(estimate),': ']);
    if isempty(MSHr)
        MSHr=estimate;
    end
    GRDNGr = input(['Give grading for faces ',num2str(iz-1),' - ',num2str(iz+ir-2),',default is [ones] (in [ ... ] format): ']);
    if isempty(GRDNGr)
        GRDNGr = ones(1,ir);
    end
    while length(GRDNGr)~=ir
        GRDNGr = input(['Icorrect array length!! Give grading for faces ',num2str(iz-1),' - ',num2str(iz+ir-2),',default is [ones] (in [ ... ] format): ']);
    end
    if ir>1
        MSHr(2)= round(MSHr(1)*(r(2)-r(1))/r(1)/GRDNGr(2)^expon)+1;
    end
    for i=3:ir
        MSHr(i)= round(MSHr(i-1)*(r(i)-r(i-1))/(r(i-1)-r(i-2))/GRDNGr(i)^expon)+1;
    end
    disp('Estimated radial mesh distribution:')
    disp(MSHr)
    temp = input('Update radial mesh distribution, or press ENTER to continue: ');
    if ~isempty(temp)
        MSHr=temp;
    end

end

function [ Bs,count ] = Boundaries(  z,r,B  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    iz=length(z);
    ir=length(r);
    BOU = blockBounds(z,r);
    figure;    
    axis([min(z)-0.1*abs(max(z)) , max(z)+0.1*abs(max(z)), 0 , max(r)+0.1*abs(max(r))]);
    hold;
    
    count = 1;
    
% Left Boundaries
    countL = 1;
    L=[];
    for j=0:ir-1
        base=j*(iz-1);
        i=0;
        if isempty(find(B==base+i));
            if base==0
                [L,countL] = base0L(L,countL,i);
            else
                [L,countL] = base1L(L,countL,base,i,iz);
            end
            plot([BOU(base+i+1,10),BOU(base+i+1,10)],BOU(base+i+1,11:12),'color','black');
            text(BOU(base+i+1,10),(BOU(base+i+1,11)+BOU(base+i+1,12))/2,num2str(count),'Color','red');
            count = count+1;
        end
        for i=1:iz-2
            if ~isempty(find(B==base+i-1)) && isempty(find(B==base+i));
                if base==0
                    [L,countL] = base0L(L,countL,i);
                else
                    [L,countL] = base1L(L,countL,base,i,iz);
                end
                plot([BOU(base+i+1,10),BOU(base+i+1,10)],BOU(base+i+1,11:12),'color','black');
                text(BOU(base+i+1,10),(BOU(base+i+1,11)+BOU(base+i+1,12))/2,num2str(count),'Color','red');
                count = count+1;
            end
        end
    end
    
% Right Boundaries
    countR = 1;
    R=[];
    for j=0:ir-1
        base=j*(iz-1);
        for i=1:iz-2
            if isempty(find(B==base+i-1)) && ~isempty(find(B==base+i));
                if base==0
                    [R,countR] = base0R(R,countR,i);
                else
                    [R,countR] = base1R(R,countR,base,i,iz);
                end
                plot([BOU(base+i,4),BOU(base+i,4)],BOU(base+i,5:6));
                text(BOU(base+i,4),(BOU(base+i,5)+BOU(base+i,6))/2,num2str(count),'Color','red');
                count = count+1;
            end
        end
        i = iz-1;
        if isempty(find(B==base+i-1));
            if base==0
                [R,countR] = base0R(R,countR,i);
            else
                [R,countR] = base1R(R,countR,base,i,iz);
            end
            plot([BOU(base+i,4),BOU(base+i,4)],BOU(base+i,5:6));
            text(BOU(base+i,4),(BOU(base+i,5)+BOU(base+i,6))/2,num2str(count),'Color','red');
            count = count+1;
        end
    end
    
% Up Boundaries
    countU = 1;
    U=[];
    for j=1:iz-1
        base=(ir-1)*(iz-1)+j-1;
        i=0;
        if isempty(find(B==base-i));
            [U,countU] = base1U(U,countU,ir-i,j,iz);
            plot(BOU(base+1,7:8),[BOU(base+1,9),BOU(base+1,9)],'color','black');
            text((BOU(base+1,7)+BOU(base+1,8))/2,BOU(base+1,9),num2str(count),'Color','red');
            count = count+1;
        end
         for i=1:ir-1
             if ~isempty(find(B==base-(i-1)*(iz-1))) && isempty(find(B==base-i*(iz-1)));
                [U,countU] = base1U(U,countU,ir-i,j,iz);
                plot(BOU(base-i*(iz-1)+1,7:8),[BOU(base-i*(iz-1)+1,9),BOU(base-i*(iz-1)+1,9)],'color','black');
                text((BOU(base-i*(iz-1)+1,7)+BOU(base-i*(iz-1)+1,8))/2,BOU(base-i*(iz-1)+1,9),num2str(count),'Color','red');
                count = count+1;
             end
         end
     end
    
% Down Boundaries
    countD = 1;
    D=[];
    for j=1:iz-1
        base=(ir-1)*(iz-1)+j-1;
        for i=1:ir-1
            if isempty(find(B==base-(i-1)*(iz-1))) && ~isempty(find(B==base-i*(iz-1)));
               [D,countD] = base1D(D,countD,ir-i,j,iz);
               plot(BOU(base-(i-1)*(iz-1)+1,1:2),[BOU(base-(i-1)*(iz-1)+1,3),BOU(base-(i-1)*(iz-1)+1,3)],'color','red');
               text((BOU(base-(i-1)*(iz-1)+1,1)+BOU(base-(i-1)*(iz-1)+1,2))/2,BOU(base-(i-1)*(iz-1)+1,3),num2str(count),'Color','red');
               count = count+1;
            end
        end
    end
    
    Bs = [L;R;U;D];
    count = count-1;
%L
%R
%U
%D
end

function [bou,count] = base0L(bou,count,i)
    bou(count,1:4)=i*4+3:-1:i*4;
    count = count+1;
end

function [bou,count] = base0R(bou,count,i)
    bou(count,1:4)=i*4:i*4+3;
    count = count+1;
end

function [bou,count] = base1L(bou,count,base,i,iz)
    base2=(base/(iz-1)-1)*iz*4+i*4;
    base1=base/(iz-1)*iz*4+i*4;
    bou(count,1:4)=[base1+1,base1  ,base2  ,base2+1];
    bou(count,5:8)=[base1+2,base1+1,base2+1,base2+2];
    bou(count,9:12)=[base1+3,base1+2,base2+2,base2+3];
    bou(count,13:16)=[base1  ,base1+3,base2+3,base2  ];
    count = count+1;
end


function [bou,count] = base1R(bou,count,base,i,iz)
    base1=(base/(iz-1)-1)*iz*4+i*4;
    base2=base/(iz-1)*iz*4+i*4;
    bou(count,1:4)=[base1+1,base1  ,base2  ,base2+1];
    bou(count,5:8)=[base1+2,base1+1,base2+1,base2+2];
    bou(count,9:12)=[base1+3,base1+2,base2+2,base2+3];
    bou(count,13:16)=[base1  ,base1+3,base2+3,base2  ];
    count = count+1;
end


function [bou,count] = base1U(bou,count,i,j,iz)
    base1=(i-1)*4*iz+4*(j);
    base2=(i-1)*4*iz+4*(j-1);
    bou(count,1:4)=[base1+1,base1  ,base2  ,base2+1];
    bou(count,5:8)=[base1+2,base1+1,base2+1,base2+2];
    bou(count,9:12)=[base1+3,base1+2,base2+2,base2+3];
    bou(count,13:16)=[base1  ,base1+3,base2+3,base2  ];
    count = count+1;
end
function [bou,count] = base1D(bou,count,i,j,iz)
    base1=(i-1)*4*iz+4*(j-1);
    base2=(i-1)*4*iz+4*(j);
    bou(count,1:4)=[base1+1,base1  ,base2  ,base2+1];
    bou(count,5:8)=[base1+2,base1+1,base2+1,base2+2];
    bou(count,9:12)=[base1+3,base1+2,base2+2,base2+3];
    bou(count,13:16)=[base1  ,base1+3,base2+3,base2  ];
    count = count+1;
end

function [ BOU ] = blockBounds( z,r )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    iz=length(z);
    ir=length(r);
    close all;
    
    ib = (iz-1)*ir;
    count = 1;
    for i = 1:iz-1
        BOU(count,1:3) = [z(i),z(i+1),0];
        BOU(count,4:6) = [z(i+1),0,r(1)];
        BOU(count,7:9) = [z(i),z(i+1),r(1)];
        BOU(count,10:12) = [z(i),0,r(1)];
        count = count + 1; 
    end
    
    for i = 1:ir-1
        for j = 1:iz-1
            BOU(count,1:3) = [z(j),z(j+1),r(i)];
            BOU(count,4:6) = [z(j+1),r(i),r(i+1)];
            BOU(count,7:9) = [z(j),z(j+1),r(i+1)];
            BOU(count,10:12) = [z(j),r(i),r(i+1)];
            count = count + 1;
        end
    end
end


