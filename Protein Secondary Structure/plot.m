% Protein Secondary Structure

%% SEQUENCE AB063105
rna_seq = 'GCGGCAACAGCGGGGCCGAUGUGUAGUUGGUGACUGCCUCUCCAGAUGCUGAGGUGCCUGUAUCAUUGGCACAGGCCAGUGCUGAACCGUAGGUGGAGUAGGCUGUGCCUUCUGAAGCAGUAUCUAUUCACAAUGAAGUUGCAGUCUCCCGAAUUCCAGUCACUUUUCACAGAAGGACUGAAGAGUCUGACAGAAUUAUUUGUCAAAGAGAAUCACGAAUUAAGAAUAGCAGGAGGAGCAGUGAGGGAUUUAUUAAAUGGAGUAAAGCCUCAGGAUAUAGAUUUUGCCACCACUGCUACCCCUACUCAAAUGAAGGAGAUGUUUCAGUCGGCUGGGAUUCGGAUGAUAAACAACAGAGGAGAAAAGCACGGAACAAUUACUGCCAGGCUUCAUGAAGAAAAUUUUGAGAUUACUACACUACGGAUUGAUGUCACCACUGAUGGAAGACAUGCUGAGGUAGAAUUUACAACUGACUGGCAGAAAGAUGCGGAACGCAGAGAUCUCACUAUAAAUUCUAUGUUUUUAGGUUUUGAUGGCACUUUAUUUGACUACUUUAAUGGUUAUGAAGAUUUAAAAAAUAAGAAAGUUAGAUUUGUUGGACAUGCUAAACAGAGAAUACAAGAGGAUUAUCUUAGAAUUUUAAGAUACUUCAGGUUUUAUGGGAGAAUUGUAGACAAACCUGGUGACCAUGAUCCUGAGACUUUGGAAGCAAUUGCAGAAAAUGCAAAAGGCUUGGCUGGAAUAUCAGGAGAAAGGAUUUGGGUGGAACUGAAAAAAAUUCUUGUUGGUAACCAUGUAAAUCAUUUGAUUCACCUUAUCUAUGAUCUUGAUGUGGCUCCUUAUAUAGGUUUACCUGCUAAUGCAAGUUUAGAAGAAUUUGACAAAGUCAGUAAAAAUGUUGAUGGUUUUUCACCAAAGCCAGUGACUCUUUUGGCCUCAUUAUUCAAAGUACAAGAUGAUGUCACAAAAUUGGAUUUGAGGUUGAAGAUCGCGAAAGAGGAGAAAAACCUUGGCUUAUUUAUAGUUAAAAAUAGGAAAGAUUUAAUUAAAGCAACAGAUAGUUCAGACCCAUUGAAACCCUAUCAAGACUUCAUUAUAGAUUCUAGGGAACCUGAUGCAACUACUCGUGUAUGUGAACUACUGAAGUACCAAGGAGAGCACUGUCUCCUAAAGGAAAUGCAGCAGUGGUCCAUUCCUCCAUUUCCUGUAAGUGGCCAUGACAUCAGAAAAGUGGGCAUUUCUUCAGGAAAAGAAAUUGGGGCUCUAUUACAACAGUUGCGAGAACAGUGGAAAAAAAGUGGUUACCAAAUGGAAAAAGAUGAACUUCUGAGUUACAUAAAGAAGACCUAAAACUGAUGGCUACUAAAAAGCAGAGCAUUU';
rna_str = '....(..(...))..(....)(.(((((.(...)(..((((((((.(((((((((..((((..(..((...))).....(..(((((....)((.....)..)).)..)((((((.......))..((((((..((...(..).)))((..(..(((((...))(((((.....(..))))).))(((...))(..)..(((((...))....))..(..)))))..(..(((.((.((...)((((..(..(..)..(((...((..((....(..)))..))).))))).))....))))..))(...((((..(....)))))))....))))..)..)..(.(..))))...).)).)))).))..).))(..).))..)))(((((..))))(((..)))((.....(..(((....)..))..))).))))...))))((.(..))))))(..((..)((((((...((...))(..(((..(((.......(..))).)(..))..))..))((((..(.(((...(...)((((..(...(..(((((....)((...(...))))))((..(.((..(...)))).)))))(((....)((((.((..((.(((......)))..)(..))((((...)))(((((((..(...)(...(....)))))))(((...)))))))))))(...)))).))(((..))))))))...)))))...)))..)(((((..((..((((((.((((..((...))...)))(..))))).((.(...)(..)..))))..))))))).)(((..((..(..)((.((((...((((......((((..(((((....)).))((...))))((((...))).(.....)(...((((..((((((((.((((.(((.(...)(((((((..)(((((.....))(..)((.(...).)(.((((..(..))))).).)))))(....).))))))).(.....)))))..(..(((.......)))..))))))((((((..))))..(.(.....))))).)))).....)))((..((.(...)((......(..((..)((..))((...))))..((..((..(((..)).....).(..)))))))))..(((.(((((..((((((.((...(...)...))))))(((((((.(((..(.((..(..).)))))).))).)))....((((((..)..))))))))))...(((....(((...))..)))))(...)).....))))))))))))...)))(((((.(..))))(...)))))))))))))))))..)))))...(....)..))((....(...)))';

rnaplot(rna_str, 'Sequence', rna_seq , 'Format', 'Diagram', 'colorby', 'pair');
title('Sequence AB063105: Secondary Structure - Diagram')

rnaplot(rna_str, 'Sequence', rna_seq, 'Format', 'Graph', 'colorby', 'pair');
title('Sequence AB063105: Secondary Structure - Graph')

%% SEQUENCE tRNA (Homo sapiens tRNA-Ala (anticodon AGC) 3-1 (TRA-AGC3-1))

rna_seq = 'GGGGGUGUAGCUCAGUGGUAGAGCGCGUGCUUAGCAUGUACGAGGUCCCGGGUUCAAUCCCCGGCACCUCCA';
rna_str_esp = '(((((((..((((.......)))).(((((.......))))).....(((((.......)))))))))))).';
rna_str_obt = '..(((((.(..)))((((..((.....(..(..((..(..)((..))(..))..))))))(..))))).)).';

rnaplot(rna_str_esp, 'Sequence', rna_seq, 'Format', 'Diagram', 'colorby', 'pair');
title('tRNA Sequence: Expected Secondary Structure - Diagram')

rnaplot(rna_str_obt, 'Sequence', rna_seq, 'Format', 'Diagram', 'colorby', 'pair');
title('tRNA Sequence: Obtained Secondary Structure - Diagram')

%% ARTIFICIAL SEQUENCE

rna_seq = 'CCCCCCCCUUAAAAAAAAGCGUUUUUUUUCCGGGGGGGG';
rna_str = '((((((((..((((((((...))))))))..))))))))';

rnaplot(rna_str, 'Sequence', rna_seq , 'Format', 'Diagram', 'colorby', 'pair');
title('Artificial Sequence: Secondary Structure - Diagram')
