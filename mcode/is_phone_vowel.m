function r=is_phone_vowel(ph)
vowels={'aa','ae','ah','ao','ax','aw','ay','ea','eh','er','ey','iy','ih','oh','ow','oy','ua','uh','uw'};
r=~isempty(fsic(vowels,lower(ph)));
return