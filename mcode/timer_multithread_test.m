function timer_multithread_test()
tm = timer('StartFcn', @hello, 'TimerFcn', @hello);
start(tm);


for i1 = 1 : 3
    fprintf(1, 'Lalala\n');
end
waitfor(tm);
return

