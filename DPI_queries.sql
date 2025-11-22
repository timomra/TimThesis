# -- Find all the compounds whos name is like diphenyliodonium

SELECT 
    rec.compound_name
FROM chembl_id_lookup AS lookup
JOIN compound_records AS rec
	ON lookup.entity_id = rec.molregno
WHERE rec.compound_name LIKE "%diphenyliodonium%"
;




# -- The CHEMBL ID for DPI (and various salt forms of this compound) were found on the website. 
# -- Finding the UNIQUE molregnos assoc with DPI -- # 
SELECT 
    DISTINCT rec.molregno
    #DISTINCT lookup.chembl_id
FROM chembl_id_lookup AS lookup
JOIN compound_records AS rec
	ON lookup.entity_id = rec.molregno
WHERE lookup.chembl_id IN ("CHEMBL365739","CHEMBL397686","CHEMBL4167328")
;


-- SELECTING activity data for those molregnos 
SELECT
	 act.pchembl_value,act.standard_type,act.standard_units,
     act.standard_value,
     act.activity_comment,
     assays.description,
     tar.pref_name
FROM activities as act
JOIN assays
	ON act.assay_id = assays.assay_id
JOIN target_dictionary AS tar
	ON tar.tid=assays.tid 
WHERE 
	act.molregno IN (322987, 391098, 2270359)
;


# -- select the cannonical smiles of DPI 
SELECT * 
FROM compound_structures
WHERE molregno LIKE 322987
;

