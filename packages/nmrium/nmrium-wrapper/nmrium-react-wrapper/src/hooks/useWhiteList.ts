import { useEffect, useState, useTransition } from 'react';

async function readAllowedOrigins() {
  return fetch(
    'https://raw.githubusercontent.com/NFDI4Chem/nmrium-react-wrapper/main/src/allowed-origins.json',
  ).then((response) => response.json());
}

export function useWhiteList() {
  const [allowedOrigins, setAllowedOrigins] = useState([]);
  const [isFetchAllowedOriginsPending, startTransition] = useTransition();

  useEffect(() => {
    startTransition(() => {
      void readAllowedOrigins().then((whitelist) =>
        setAllowedOrigins(whitelist),
      );
    });
  }, []);

  return {
    allowedOrigins,
    isFetchAllowedOriginsPending,
  };
}
