function res = EF_num2Sub(inp)

%converts a numeric index to a name index and vice versa

subs =   {'ef_071411' 'ef_072111' 'ef_083111' 'ef_091211' 'ef_091311' ...
    'ef_091511' 'ef_091911' 'ef_092111' 'ef_092211' 'ef_092711' ...
    'ef_092911'  'ef_100511' 'ef_101411' 'ef_032912' 'ef_040412' ...
    'ef_040512' 'ef_040712' 'ef_040712_2' 'ef_041112' 'ef_042912' 'ef_050112'};


if isnumeric(inp)
    res = subs{inp};
elseif ischar(inp)
    res = find(strcmp(inp, subs));
else
    error ('input must be an integer or char')
end