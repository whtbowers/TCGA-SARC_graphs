library(RSQLite)

#dbConnect creates new database

# mydb <- dbConnect(RSQLite::SQLite(), "my-db.sqlite")
# dbdisconnect(mydb)
# unlink("my-db.sqlite")

# If just temp , "" for on-disk db, :memory:, or file::memory: for an in-memory database
# A temporary database will be erased as soon as it's disconnected
# mydb <- dbConnect(RSQLite::SQLite(), "")
# dbDisconnect(mydb)

mydb <- dbConnect(RSQLite::SQLite(), "")
dbWriteTable(mydb, "mtcars", mtcars)
dbWriteTable(mydb, "iris", iris)
dbListTables(mydb)

# Queries - Standard SQL 

dbGetQuery(mydb, 'SELECT * FROM mtcars LIMIT 5')

# Some variable name don't fit standard SQL structure and may need to be escaped

dbGetQuery(mydb, 'SELECT * FROM iris WHERE "Sepal.Length" < :x',
           params = list(x = 4.6))

# dbFetch fetches results in batches which don't fit in the memory

rs <- dbSendQuery(mydb, 'SELECT * FROM mtcars')
while (!dbHasCompleted(rs)) {
  df <- dbFetch(rs, n =10)
  print(nrow(df))
}

## Multiple parameterised queries
rs <- dbSendQuery(mydb, 'SELECT * FROM iris WHERE "Sepal.Length" < :x')
dbBind(rs, param = list(x = 4.5))
nrow(dbFetch(rs))

dbBind(rs, param = list(x = 4))
nrow(dbFetch(rs))

dbClearResult(rs)

# Pass multiple parameters in one call to dbBind
rs <- dbSendQuery(mtdb, 'SELECT * FROM iris WHERE "Sepal.Length = :x')
dbBind(rs, param = list(x = seq(4, 4.4, by = 0.1)))
nrow(dbFetch(rs))
dbClearResults(rs)

## Statements
dbExecute(mydb, 'DELETE FROM iris WHERE "Sepal.Length" < 4')

rs <- dbSendStatement(mydb, 'DELETE FROM iris WHERE "Sepal.Length" < :x')
dbBind(rs, param = list(x = 4.5))
dbGetRowsAffected(rs)
dbClearResult(rs)
