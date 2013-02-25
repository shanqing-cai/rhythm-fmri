function ivc = is_vowel_cmu(phn)
vowels = {'AA', 'AE', 'AH', 'AO', 'AW', 'AY', 'EH', 'ER', 'EY', 'IH', 'IY', 'OY', 'OW', 'UH', 'UW'};
      
ivc = ~isempty(fsic(vowels, phn));
return