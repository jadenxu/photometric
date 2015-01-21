clear all;
close all;


path(path, './Includes');

findex = '../build/final_normal_index.txt'; fname = '../build/final_normal.txt';


file = load(fname);
M = load(findex);
rows = M(1);
cols = M(2);

image = zeros(rows, cols, 3);

for x = 1:rows
    for y = 1:cols
        image(x, y, 1:3) = file((x-1)*cols + y, :);
    end
end

%[slant, tilt] = grad2slanttilt(image(:,:,1), image(:,:,2));
%slant=image(:,:,1)./(image(:,:,3)+0.000001); tilt=image(:,:,2)./(image(:,:,3)+0.0000001);
slant=image(:,:,1); tilt=image(:,:,2);
z3=shapeletsurf(slant,tilt,6,0.4,2,'slanttilt');
figure;
surf(z3)


