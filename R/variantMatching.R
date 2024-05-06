variant_match <- function(effect_allele, other_allele, ALT, REF, VT) {

  # Palindromic variants:
  # A-T
  # vs A-T : match (1)
  # vs T-A : flip (2)
  # anything else : wrong (3)

  # non-Palindromic variants:
  # A-G
  # vs A-G : match (1)
  # vs G-A: flip (2)
  # vs C-T: flip (5)
  # vs T-C : switch (4)
  # anything else : wrong (6)

  # Not Found == 9L
  # Match == 1L
  # Flip == 2L
  # Switch == 3L
  # Wrong == 4L


  other_allele_switched <- switch_alleles_vectorized(other_allele)
  effect_allele_switched <- switch_alleles_vectorized(effect_allele)



  if (VT == 1) {
    ##################
    ## SNP variants ##
    ##################

    palindromic_variant <- ifelse(effect_allele == other_allele_switched,
                                  TRUE,
                                  FALSE
    )

    # Not found in database
    if (is.na(REF) | REF == "")
      return(list(9L, palindromic_variant))


    if (effect_allele == ALT && other_allele == REF)
																	   
      return(list(1L, palindromic_variant))


    if(palindromic_variant)
    {
      if (effect_allele == REF && other_allele == ALT)
																		  
        return(list(2L, palindromic_variant))
      else
        return(list(4L, palindromic_variant))
    }
    else
    {
      if (other_allele_switched == REF && effect_allele_switched == ALT)
																							
        return(list(3L, FALSE))
      else if ((effect_allele == REF && other_allele == ALT) |
               (effect_allele_switched == REF && other_allele_switched == ALT))
        return(list(2L, FALSE))
      else
        return(list(4L, FALSE))
    }
  }

  else

  {
    #####################
    # nonSNP variants  ##
    #####################

    # Not found in database
    if (is.na(REF) | REF == "")
      return(list(9L, FALSE))


    if (effect_allele == "R" |
        (effect_allele == REF && other_allele == ALT) |
        (effect_allele_switched == REF && other_allele_switched == ALT))
      return(list(2L, FALSE))


    if (other_allele_switched == REF && effect_allele_switched == ALT) {
																						   
      return(list(3L, FALSE))
    }

    if(effect_allele %in% c("R","D","I") | (effect_allele == ALT && other_allele == REF))
																										 
      return(list(1L, FALSE))


    return(list(4L, FALSE))

  }

}
