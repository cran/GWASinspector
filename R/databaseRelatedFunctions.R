#' Check the data in the reference database
#'
#' This function displays the summary of the database, including how many tables are in the database file, number of data rows for each data table and the first row of each table.
#'
#' @param path character. full path to the database file (*.sqlite)
#' @return This function returns a data table including the summary of the specified database. This is neccessary to check the consistency and validity of an unknown or new database file.
#' @examples
#' check.database(system.file("extdata", "sample_db.sqlite", package = "GWASinspector"))
#' @note
#' First column include the names of the tables in database
#'
#' Second column is the number of rows in each table
#'
#' Next columns are the first row of each table
#'
check.database <- function(path)
{

  if(tools::file_ext(path) != 'sqlite')
    stop('File is not sqlite type database!')

  if (!file.exists(path))
    stop('File not found!')

  # message(paste('Reading database file:','"' ,basename(path),'"'))

  check<-tryCatch(
    {
      mydb<-RSQLite::dbConnect(RSQLite::SQLite(),path)
      tables<-RSQLite::dbListTables(mydb)
      # pb <- txtProgressBar(min = 0, max = length(tables), style = 3,width = 65)
      # counter <- 1


      ### ==========================
      rowcount <- data.table()
      rowTable <- data.table()

      for(i in tables){
        query <- sprintf('select count(*) from %s',i)
        rowcount <- rbind(rowcount,
                          as.data.table(RSQLite::dbGetQuery(mydb, query)))

        rowTable <- rbind(rowTable,
                          as.data.table(RSQLite::dbGetQuery(mydb, sprintf('select * from %s LIMIT 1',i))), fill = TRUE)

        # setTxtProgressBar(pb, counter)
        # counter <- counter + 1
      }

      colnames(rowcount) <- 'number of rows'
      RSQLite::dbDisconnect(mydb)
      #close(pb)

      ### ==========================


      tables <- as.data.table(tables)
      colnames(tables) <- 'table name'
      ### ==========================
      number <- as.data.table(1:nrow(tables))
      colnames(number) <- '#'

      ### ==========================
      # message('\n- first column is tables of the database')
      # message( '- second column is the number of rows in that table')
      # message( '- next columns are the first row of each table')

      tbl <- cbind(number,tables,rowcount,rowTable)
     # print(knitr::kable(as.data.table(tbl),format = 'rst',align = 'l'))
      return(as.data.table(tbl))

    },error=function(err)
    {
      stop(paste('Error reading database file:',err$message))
    })

}
