function IntMap = IntensityPlaneAnyDrop( rayData, imgDim, res )
%This is for the parallelized collection of intensity, with the hope of
%being able to collect angular data, intensity data all in parallel for a
%fully characterized single experiment. 

%Here in multipartrace, we grab the rays final data, location, direction,
%and amplitude. These are then used in order to place them over a matrix
%below with improfile, which then rounds them in order to estimate their
%path and thus 'intensity'. 

%The main goal is to reduce either memory writing of ray cells for large
%ray simulations, and also to reduce expensive computation of other methods

%rayData is a n by 5 array where (x loc) (y loc) (x dir) (y dir) (amplitude)
%This is specialized to handle only last ray data, where we collect it
%directly after the droplet has been traced, and therefore can just apply
%these rays over the map.

xs = imgDim(1)/(2);
ys = imgDim(2)/(2);
img = zeros(imgDim);
img2 = img;
ext = imgDim(1);
for i = 1:size(rayData,1)
    x1 = rayData(i,1)*res + xs;
    x2 = rayData(i,1)*res + rayData(i,3)*ext + xs;
    y1 = rayData(i,2)*res + ys;
    y2 = rayData(i,2)*res + rayData(i,4)*ext + ys;
    a = rayData(i,5);
    %
    %d = round(sqrt((x1-x2)^2+(y1-y2)^2));
    %cx = linspace(x1,x2,d)';
    %cy = linspace(y1,y2,d)';
    %
    %Whats written between the two % above replaces improfile because its
    %more resource efficient and should get the same job done
    [cx,cy,~] = improfile(img,[x1 x2],[y1 y2],'bicubic');
    line = round([cx cy]);
    data = unique(line,'rows');
    data = data(data(:,1) > 0, :);
    data = data(data(:,2) > 0, :);
    data = data(data(:,2) <= imgDim(1), :);
    data = data(data(:,1) <= imgDim(2), :);
    
        for k = 1:(size(data,1)) %for k = 1:(size(data,1) -1)
            img2(data(k,2),data(k,1)) = img2(data(k,2),data(k,1)) + a;
        end

end
    IntMap = img2;       
end

