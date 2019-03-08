function [blurlen, bluranle] = auto2motion(text_aut)
ifshowfigure = 0;
text_aut_max = text_aut;
text_aut_max(text_aut_max==1) = 0;
% text_aut_max(text_aut_max~=1) = 0;
text_aut_max(text_aut_max~=0) = 1;
p=FastPeakFind(text_aut_max,0);
if ifshowfigure==1
figure,imshow(text_aut_max); hold on

plot(p(1:2:end),p(2:2:end),'r+','MarkerSize',15);%title('find peak points and lines');
end
point_peak1 = [p(1:2:end) p(2:2:end)];

text_aut_peak = text_aut_max;
if size(point_peak1,1)==2
    endpoints = point_peak1;
    mm = endpoints(2,:)-endpoints(1,:);
    % xmm = norm(mm(:));
    bluranle = -atand(mm(2)/mm(1));
    blurlen = ceil(norm(mm)/2);
    if ifshowfigure==1
    figure, imshow(text_aut_max);hold on
    plot(endpoints(:,1),endpoints(:,2),'LineWidth',4,'Color','r');
    set(gca,'FontWeight','bold','FontSize',20)
    saveas(gca,('./result/path.epsc'))
    end
else
    % line
    [H, theta, rho]= hough(text_aut_peak,'RhoResolution', 1);
    peak=houghpeaks(H,5);
    lines=houghlines(text_aut_peak,theta,rho,peak);
    
    TF = isempty(lines);
    
    % figure,imshow(text_aut_max); hold on
    if TF ~=1
        meanxy = zeros(2,2);
        meanangl = 0;
        %         figure,imshow(text_aut_peak); hold on
        for layer=1:length(lines)
            xy=[lines(layer).point1;lines(layer).point2];
            if ifshowfigure==1
            plot(xy(:,1),xy(:,2),'LineWidth',4,'Color','r');
            end
            meanxy = meanxy + xy;
            meanangl = meanangl+lines(layer).theta;
            
            %         figure,imshow(text_aut_peak,[]),title('move line'),hold on
        end
        % end
        if ifshowfigure==1
        set(gca,'FontWeight','bold','FontSize',20)
        saveas(gca,('./result/peakline.epsc'))
        end
        
        meanxy = floor(meanxy./layer);
        meanangl = floor(meanangl./layer);
        
        
        rankangl = zeros(layer,1);
        for layer=1:length(lines)
            
            rankangl(layer) = abs(meanangl - lines(layer).theta);
        end
        [~,idxmin] = min(rankangl(:));
        
        endpoints = [lines(idxmin).point1;lines(idxmin).point2];
        estangle = lines(idxmin).theta;
        if meanangl>0
            bluranle = 90-abs(estangle);
        else
            bluranle = -(90-abs(estangle));
        end
        %     endpoints = intersect(point_peak1,point_peak2,'rows');
        if ifshowfigure==1
        figure, imshow(text_aut_max);hold on
        plot(endpoints(:,1),endpoints(:,2),'LineWidth',4,'Color','r');
        set(gca,'FontWeight','bold','FontSize',20)
        saveas(gca,('./result/path.epsc'))
        end
        mm = endpoints(end,:)-endpoints(1,:);
        blurlen = ceil(norm(mm)/2);
    else
        [x,y] = find(text_aut_peak==1); 
        nx =sort(x);mx = abs(nx(1)-nx(end));ny =sort(y);my = abs(ny(1)-ny(end));
        if (mx>my)&&(mx~=0)
            [~,idxminx] = min(x);
            [~,idxmaxx] = max(x);
            endpoints = [y(idxminx) min(x); y(idxmaxx) max(x)];
            if ifshowfigure==1
            figure, imshow(text_aut_max);hold on
            plot(endpoints(:,1),endpoints(:,2),'LineWidth',4,'Color','r');
            set(gca,'FontWeight','bold','FontSize',20)
            saveas(gca,('./result/path.epsc'))
            end
            mm = endpoints(2,:)-endpoints(1,:);
        else
            [~,idxminx] = min(y);
            [~,idxmaxx] = max(y);
            endpoints = [ min(y) x(idxminx) ;y(idxmaxx) x(idxmaxx)];
            if ifshowfigure==1
            figure, imshow(text_aut_max);hold on
            plot(endpoints(:,1),endpoints(:,2),'LineWidth',4,'Color','r');
            set(gca,'FontWeight','bold','FontSize',20)
            saveas(gca,('./result/path.epsc'))
            end
            mm = endpoints(2,:)-endpoints(1,:);
        end
        
        % xmm = norm(mm(:));
        bluranle = -atand(mm(2)/mm(1));
        blurlen = ceil(norm(mm)/2);
    end
end