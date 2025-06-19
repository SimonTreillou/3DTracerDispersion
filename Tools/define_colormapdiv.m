function [cmap] = define_colormapdiv(c1,c2,cols,perc)
    % Define colormap
    % c1: positive color
    % c2: negative color
    % cols: number of columns
    % perc: position of white
    c1=c1(1:3); c2=c2(1:3);
    col1=int16((1-perc)*cols);
    col2=int16((perc)*cols);


    grad=flip(linspace(0,1,col1));
    C1=repmat(c1,[col1,1]);
    C1=C1 + (([1,1,1]-C1).*grad');

    grad=linspace(0,1,col2);
    C2=repmat(c2,[col2,1]);
    C2=C2 + (([1,1,1]-C2).*grad');

    cmap=cat(1,C2,C1);
end