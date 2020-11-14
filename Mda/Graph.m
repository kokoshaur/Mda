for i = 0:9
    spectrum = load(strcat('D:\1)Me\Ўарага\5 семак\цос\Mda Ч копи€ (2)\Mda\after\MDA' , num2str(i),'.txt'));
    x = spectrum(:,1);
    spectrum2 = load(strcat('D:\1)Me\Ўарага\5 семак\цос\Mda Ч копи€ (2)\Mda\before\ansver' , num2str(i),'.txt'));
    y = spectrum2(:,1);
    
    figure(i+1);
    subplot(1,2,1); plot(y,'b');
    title(strcat('Before #', num2str(i+1)));
    
    subplot(1,2,2); plot(x,'r');
    title(strcat('After #', num2str(i+1)));
end