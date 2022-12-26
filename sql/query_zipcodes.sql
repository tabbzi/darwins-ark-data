SELECT user_id,
       MAX(CASE WHEN meta_key = 'billing_state' THEN meta_value ELSE NULL END) AS billing_state,
       MAX(CASE WHEN meta_key = 'billing_country' THEN meta_value ELSE NULL END) AS billing_country,
       MAX(CASE WHEN meta_key = 'billing_postcode' THEN meta_value ELSE NULL END) AS billing_postcode,
       MAX(CASE WHEN meta_key = 'state' THEN meta_value ELSE NULL END) AS state,
       MAX(CASE WHEN meta_key = 'country' THEN meta_value ELSE NULL END) AS country,
       MAX(CASE WHEN meta_key = 'postcode' THEN meta_value ELSE NULL END) AS postcode,
       MAX(CASE WHEN meta_key = 'shipping_state' THEN meta_value ELSE NULL END) AS shipping_state,
       MAX(CASE WHEN meta_key = 'shipping_country' THEN meta_value ELSE NULL END) AS shipping_country,
       MAX(CASE WHEN meta_key = 'shipping_postcode' THEN meta_value ELSE NULL END) AS shipping_postcode
FROM `wp_usermeta`
WHERE meta_key IN ('billing_state','billing_country','billing_postcode','state','country','postcode','shipping_state','shipping_country','shipping_postcode')
GROUP BY user_id;

