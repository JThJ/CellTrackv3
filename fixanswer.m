function [answer] = fixanswer(answer)

% if char contains numbers this conversion is needed
[h,w]=size(answer);
if h>w
    answer = answer';
end