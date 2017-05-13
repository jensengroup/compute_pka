N_sp3_smiles = ( ('[NH3+]c1ccccc1',4.6,'aniline_1'),
               ('C[NH2+]c1ccccc1',4.9,'ethylaniline_1'),
               ('C[NH+](C)c1ccccc1',5.9,'diethylaniline_1'),
               ('C[NH3+]',10.6,'ethylamine_1'),
               ('C[NH2+]C',11.0,'diethylamine_1'),
               ('C[NH+](C)C',10.9,'triethylamine_1'),
               ('CNS(C)(=O)=O',13.3,'ethylmethanesulfonamide_1') )

#should add benzamidine
N_sp2_smiles = ( ('C=[NH2+]', 10.9,'ethanimine_1'),
                 ('C=[NH+]C', 11.0,'ethylethanimine_1'),
                 ('CC(N)=[NH2+]', 13.0,'propionimidamide_1'),
                 ('NC(=[NH2+])', 12.9,'ethylpropionimidamide_1'),
                 ('C[NH+]=C(N)C', 12.9,'ethylpropionimidamide_1'),
                 ('CNC(CC)=[NH+]C', 11.9,'diethylpropionimidamide_1'),
                 ('CC(=[NH2+])Nc1ccccc1', 10.3,'phenylethanimidamide_1'),
                 ('CNC(N)=[NH2+]', 12.8,'ethylguanidine_1'),
                 ('CNC(=[NH2+])NC', 13.0,'diethylguanidine_1'),
                 ('CNC(N)=[NH+]C', 13.0,'diethylguanidine_1'),
                 ('CNC(NCC)=[NH+]C', 12.8,'triethylguanidine_1'))

N_heterocycle_smiles = ( ('c1c[nH+]c[nH]1', 7.0,'imidazole_1'),
                         ('c1ccc2[nH+]c[nH]c2c1', 5.6,'benzimidazole_1'),
                         ('c1coc[nH+]1', 0.8,'oxazole_1'),
                         ('c1csc[nH+]1', 2.5,'thiazole_1'),
                         ('c1ccc2sc[nH+]c2c1', 0.8,'benzothiazole_1'),
                         ('c1c[nH+]ccn1', 0.4,'pyrazine_1'),
                         ('c1c[nH][nH+]c1', 2.8,'pyrazole_1'),
                         ('c1cc[nH+]nc1', 2.5,'pyridazine_1'),
                         ('c1cc[nH+]cc1', 5.2,'pyridine_1'),
                         ('c1ccc2[nH+]cccc2c1', 4.8,'quinoline_1'),
                         ('c1cnc[nH+]c1', 1.1,'pyrimidine_1'))

O_sp2_smiles = ( ('Oc1ccccc1', 10.0,'phenol_0'),
                 ('C=CO', 9.9,'ethenol_0'),
                 ('C=C(C)O', 11.6,'propen2ol_0'),
                 ('C(=O)O', 4.8,'aceticacid_0'),
                 ('O=C(O)c1ccccc1', 4.2,'benzoicacid_0'))

S_smiles = ( ('CCS', 10.6,'ethanethiol_0'),
             ('C=CS', 9.4,'ethenethiol_0'),
             ('C=C(C)S', 9.6,'propene2thiol_0'),
             ('Sc1ccccc1', 6.6,'benzenethiol_0'))

