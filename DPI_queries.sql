# -- The CHEMBL ID for DPI (and various salt forms of this compound) were found on the website. 
# -- Finding the UNIQUE molregnos assoc with DPI -- # 
SELECT 
    DISTINCT rec.molregno
    #DISTINCT lookup.chembl_id
FROM chembl_id_lookup AS lookup
JOIN compound_records AS rec
	ON lookup.entity_id = rec.molregno
WHERE lookup.chembl_id IN ("CHEMBL365739","CHEMBL397686","CHEMBL4167328")