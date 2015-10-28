-- Export some data for fun

-- collar
select * from header order by Loc_utmx, loc_utmy LIMIT 5;
-- survey 
select * from survey where HOLE_ID in (select [HOLE-ID] from header order by Loc_utmx, loc_utmy LIMIT 5); 
-- assay 
select * from assays where HOLE_ID in (select [HOLE-ID] from header order by Loc_utmx, loc_utmy LIMIT 5); 
-- litho
select * from litho where HOLE_ID in (select [HOLE-ID] from header order by Loc_utmx, loc_utmy LIMIT 5); 
-- structure
select * from assay where HOLE_ID in (select [HOLE-ID] from header order by Loc_utmx, loc_utmy LIMIT 5); 
