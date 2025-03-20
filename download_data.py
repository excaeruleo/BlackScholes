import yfinance as yf
import pandas as pd
from datetime import datetime, timedelta

# Get today's date
end_date = datetime.now()
# Calculate start date (2 years ago)
start_date = end_date - timedelta(days=2*365)

# Download AAPL stock data
aapl = yf.Ticker("AAPL")
stock_data = aapl.history(start=start_date, end=end_date)

# Get current options data instead of historical
options_data = pd.DataFrame()
try:
    # Get all available expiration dates
    expirations = aapl.options

    # Get both calls and puts for the nearest expiration
    if expirations:
        opt = aapl.option_chain(expirations[0])
        calls_data = opt.calls[['strike', 'lastPrice', 'volume', 'impliedVolatility']]
        puts_data = opt.puts[['strike', 'lastPrice', 'volume', 'impliedVolatility']]
        
        # Save current options data separately
        calls_data.to_csv('aapl_calls.txt', index=False)
        puts_data.to_csv('aapl_puts.txt', index=False)
except Exception as e:
    print(f"Error fetching options data: {e}")

# Save stock data
stock_data['Close'].to_csv('aapl_stock.txt', header=True)

print("Stock data has been saved to aapl_stock.txt")
print("Current calls data has been saved to aapl_calls.txt")
print("Current puts data has been saved to aapl_puts.txt")
